import collections
import itertools
import os

import numpy
import pysam

OUTPUT_FILE_NAME = "intron_quantifications.txt"

class IntronQuantificationStep:
    def __init__(self, geneinfo_file_path, flank_size = 1500):
        '''
        geneinfo_file_path is a path to a tab-separated file with the following fields:
        chrom  strand  txStart  txEnd  exonCount  exonStarts  exonEnds  transcriptID  geneID  geneSymbol  biotype
        '''
        self.geneinfo_file_path = geneinfo_file_path
        self.flank_size = flank_size

    def execute(self, aligned_file_path, output_directory, forward_read_is_sense):
        with pysam.AlignmentFile(aligned_file_path, "rb") as alignments:

            # Read in the gene alignment information
            chrom_lengths = dict(zip(alignments.references, alignments.lengths))
            self.load_gene_info(self.geneinfo_file_path, chrom_lengths)

            # (chrom, strand) -> {(transcriptID, intron_number) -> read_count}
            intron_read_counts = collections.defaultdict(lambda: collections.defaultdict(lambda:0))
            # chrom -> {intron -> read count}
            intergenic_read_counts = collections.defaultdict(lambda: collections.defaultdict(lambda:0))

            unpaired_reads = dict()

            # Go through all reads, and compare
            skipped_chromosomes = []
            for read in alignments.fetch(until_eof=True):
                # Use only uniquely mapped reads with both pairs mapped
                if read.is_unmapped or not read.is_proper_pair or not read.get_tag("NH") == 1:
                    continue

                try:
                    mate = unpaired_reads[read.query_name]
                except KeyError:
                    # mate not cached for processing, so cache this one
                    unpaired_reads[read.query_name] = read
                    continue

                # Read is paired to 'mate', so we process both together now
                # So remove the mate from the cache since we're done with it
                del unpaired_reads[read.query_name]

                chrom = read.reference_name
                # Check annotations are available for that chromosome
                if chrom not in self.intergenics:
                    if chrom not in skipped_chromosomes:
                        print(f"Alignment from chromosome {chrom} skipped")
                        skipped_chromosomes.append(chrom)
                    continue

                # According to the SAM file specification, this CAN fail but I don't understand why it would
                # so just throw this assert in to verify that it doesn't, at least for now
                assert read.is_reverse != mate.is_reverse

                # Figure out the fragment's strand - depends on whether the forward or reverse reads are 'sense'
                if forward_read_is_sense:
                    strand = "-" if (read.is_reverse and read.is_read1) or (not read.is_reverse and read.is_read2) else "+"
                else:
                    strand = "+" if (read.is_reverse and read.is_read1) or (not read.is_reverse and read.is_read2) else "-"

                mintron_starts, mintron_ends = self.mintrons[(chrom, strand)]
                intergenic_starts, intergenic_ends = self.intergenics[chrom]

                mintrons_touched = set()
                intergenics_touched = set()

                # check what regions our blocks touch
                for start,end in itertools.chain(read.get_blocks(), mate.get_blocks()):
                    # NOTE: pysam is working in 0-based, half-open coordinates! We are using 1-based
                    start = start+1

                    # We may intersect this mintron, depending on where its end lies
                    # or this could be -1 if we start before all mintrons
                    last_mintron_before = numpy.searchsorted(mintron_starts, start, side="right") - 1

                    # We definitely do not intersect this mintron, it starts after our end
                    first_mintron_after = numpy.searchsorted(mintron_starts, end, side="right")
                    # but all between the two, we do intersect

                    if last_mintron_before == -1:
                        # No introns start before us
                        mintrons_touched.update(range(0, first_mintron_after))
                    elif mintron_ends[last_mintron_before] >= start:
                        # We do intersect the last one to start before us
                        mintrons_touched.update(range(last_mintron_before, first_mintron_after))
                    else:
                        # We only start intersecting the first one after that
                        mintrons_touched.update(range(last_mintron_before+1, first_mintron_after))


                    # Now do the same thing as above for intergenic regions
                    last_intergenic_before = numpy.searchsorted(intergenic_starts, start, side="right") - 1
                    first_intergenic_after = numpy.searchsorted(intergenic_starts, end, side="right")
                    if last_intergenic_before == -1:
                        intergenics_touched.update(range(0, first_intergenic_after))
                    elif intergenic_ends[last_intergenic_before] >= start:
                        intergenics_touched.update(range(last_intergenic_before, first_intergenic_after))
                    else:
                        intergenics_touched.update(range(last_intergenic_before+1, first_intergenic_after))


                # Accumulate the reads
                for mintron in mintrons_touched:
                    for intron in self.mintron_annotations[(chrom,strand)][mintron]:
                        intron_read_counts[(chrom,strand)][intron] += 1

                for intergenic in intergenics_touched:
                    intergenic_read_counts[chrom][intergenic] += 1

        normalized_intron_read_counts, percentages, gene_level_counts, transcript_level_counts \
                = normalize_read_counts(intron_read_counts, self.effective_transcript_lengths)

        # Write out the results to output file
        output_file_path =  os.path.join( output_directory, OUTPUT_FILE_NAME)
        with open(output_file_path, "w") as output_file:
            output_file.write("#gene_id\ttranscript_id\tgene_reads\ttranscript_reads\tintron_percentages\n")
            # take transcripts from all chromosomes and combine them, sorting by gene id and then transcript id
            transcripts = sorted( itertools.chain(*(info for info in self.transcript_info.values())) )
            for (gene_id, transcript_id), _, _ in transcripts:
                num_introns = len(self.intron_info[gene_id, transcript_id][0])
                transcript_reads = str(transcript_level_counts[gene_id, transcript_id])
                gene_reads = str(gene_level_counts[gene_id])
                pcts = [str(percentages[(gene_id, transcript_id), i]) for i in range(num_introns)]

                output_file.write('\t'.join([gene_id,
                                             transcript_id,
                                             gene_reads,
                                             transcript_reads,
                                             ','.join(pcts),
                                             ]) + '\n')

    def load_gene_info(self, geneinfo_file_path, chrom_lengths):
        # Assumes that geneinfo_file is sorted by feature starts within a chromosome!!
        # and that the chromosomes in geneinfo_file are all present in chrom_lengths (i.e. in the bam file header)

        # Gather exon boundaries from all annotations, sorted by their start
        # Collect exons by chromosome and strand, and intergenics just by chromosome
        make_deque = lambda: collections.deque()
        make_dict = lambda: dict()
        exons = collections.defaultdict(make_deque)
        genics = collections.defaultdict(make_deque)
        all_introns = collections.defaultdict(make_deque)
        transcript_info = collections.defaultdict(make_deque)
        gene_info = collections.defaultdict(make_dict)
        with open(geneinfo_file_path) as geneinfo_file:
            for line in geneinfo_file:
                if line.startswith("#"):
                    continue
                chrom, strand, txStart, txEnd, exonCount, exonStarts, exonEnds, transcriptID, geneID, geneSymbol, biotype = line.split('\t')
                txStart, txEnd = int(txStart), int(txEnd)
                exonStarts = [int(start) for start in exonStarts.split(',')]
                exonEnds = [int(end) for end in exonEnds.split(',')]
                exons[(chrom,strand)].extend(zip(exonStarts, exonEnds)) # Track strandedness of exons
                genics[chrom].append( (txStart, txEnd) )
                transcript_info[chrom].append( ((geneID, transcriptID), txStart, txEnd) )

                intronStarts = [end + 1 for end in exonEnds[:-1]]
                intronEnds = [start - 1 for start in  exonStarts[1:]]
                all_introns[(chrom,strand)].append( ((geneID, transcriptID), (intronStarts, intronEnds)) )

                if geneID in gene_info:
                    geneStart, geneEnd = gene_info[geneID]
                    gene_info[geneID] = (min(geneStart, txStart), max(geneEnd, txEnd))
                else:
                    gene_info[geneID] = (txStart, txEnd)

        # Convert to arrays (array[0] is array of starts, array[1] array of ends)
        exons = {chrom: numpy.sort(numpy.array(exon).T) for chrom,exon in exons.items()}
        genics = {chrom: numpy.sort(numpy.array(genic).T) for chrom,genic in genics.items()}

        # Merge regions that overlap
        merged_exons = {chrom: merge_regions(exon) for chrom, exon in exons.items()}
        merged_genics = {chrom: merge_regions(genic) for chrom, genic in genics.items()}

        # Intergenics are complements of all the genics
        intergenics = {chrom: complement_regions(genic, chrom_lengths[chrom]) for chrom,genic in merged_genics.items()}

        # Shrink each intergenic to account for flanking regions
        intergenics = {chrom: shrink_regions(intergenic, self.flank_size) for chrom, intergenic in intergenics.items()}

        # Find the flanking regions of each transcript (possibly empty...)
        flanks = {chrom: flanking_regions(transcript_info[chrom], merged_genics[chrom], self.flank_size) for chrom in genics.keys()}
        # Add flanking regions as "introns" to each gene
        def add_flank(ID, starts, ends, chrom):
            left, right = flanks[chrom][ID]
            return (ID, ([left[0]] + starts + [right[0]], [left[1]] + ends + [right[1]]))
        all_introns = {(chrom,strand): [add_flank(ID, intronStarts, intronEnds, chrom)
                                            for ID, (intronStarts, intronEnds) in introns]
                                        for (chrom,strand), introns in all_introns.items()}

        # intronic regions are ones that between these merged exons that ARENT intergenic
        # so combine exons and intergenics
        # These regions are disjoint, so just need to sort them together
        intergenics_or_exons = {(chrom,strand): numpy.concatenate((exon, intergenics[chrom]), axis=1) for (chrom,strand),exon in merged_exons.items()}
        [intergenic_or_exon.sort() for intergenic_or_exon  in intergenics_or_exons.values()] # And sort them into the right order
        # And merge them to combine adjacent regions
        intergenics_or_exons = {(chrom,strand): merge_regions(intergenic_or_exon) for (chrom,strand),intergenic_or_exon in intergenics_or_exons.items()}

        # We define mintrons to be a region that is not intergenic and is not in an exon of THIS strand
        # i.e. is in at least one gene on this strand and not in any exon of any gene of this strand
        # (for min-introns, of which there might be several in a real intron since exons of other genes could break up the intron)
        # Analogously, we should define the 'merged exons' above to be 'maxons' since they're max of all touching exons...
        mintrons = {(chrom,strand): complement_regions(regions, chrom_lengths[chrom]) for (chrom,strand),regions in intergenics_or_exons.items()}

        # Annotate the mintrons by the transcriptID + introns
        # and find the effective lengths of the introns
        mintron_annotations, effective_lengths = annotate_regions(mintrons, all_introns, gene_info)

        self.merged_exons = merged_exons
        self.mintrons = mintrons
        self.merged_genics = merged_genics
        self.intergenics = intergenics
        self.mintron_annotations = mintron_annotations
        self.effective_transcript_lengths = effective_lengths
        self.transcript_info = transcript_info

        # Merge intron info across all chromosomes
        self.intron_info = dict()
        for  chrom_introns in all_introns.values():
            self.intron_info.update(dict(chrom_introns))

        print(f"Loaded info on {sum( v.shape[1] for v in intergenics.values())} intergenic regions and {sum(v.shape[1] for v in mintrons.values())} mintrons")

def normalize_read_counts(read_counts, effective_lengths):
    normalized_counts = dict()
    transcript_level_counts = collections.Counter()
    gene_level_counts = collections.Counter()
    for chrom, chrom_read_counts in read_counts.items():
        for intron, count in chrom_read_counts.items():
            (gene_id, transcript_id), intron_num = intron
            length = effective_lengths[gene_id, transcript_id]
            norm_count = count / length * 1_000 # FPK (no M) - i.e. reads per 1000 effective bases
            normalized_counts[intron] = norm_count
            transcript_level_counts[gene_id, transcript_id] += norm_count
            gene_level_counts[gene_id] += norm_count

    # Compute the percentage of the intron reads going to this intron, among all introns of a transcript
    percentages = collections.Counter()
    for intron, count in normalized_counts.items():
        transcript, intron_num = intron
        percent = count / transcript_level_counts[transcript]
        percentages[intron] = percent
    return normalized_counts, percentages, gene_level_counts, transcript_level_counts

def merge_regions(region_arrays):
    # Given region array, i.e. an numpy array such that
    # region_array[0] lists the start positions
    # region_array[1] lists the end positions
    # Sorted by start poisition
    # Merge any regions that overlap (or are adjacent)
    last_start = -1
    last_end = -1
    merged = collections.deque()
    for start,end in region_arrays.T:
        if start <= last_end + 1:
            # Overlaps with previous start...
            # ... OR is adjacent to previous exon (want to join touching exons)
            if end > last_end:
                # but extends further
                merged[-1] = (last_start, end)
                last_end = end
            else:
                pass # Skip, contained in the existing exon
        else:
            # New exon, no overlap
            merged.append((start, end))
            last_start = start
            last_end = end

    return numpy.array(merged).T

def complement_regions(region_array, total_region_length):
    # Given region array, i.e. an numpy array such that
    # region_array[0] lists the start positions
    # region_array[1] lists the end positions
    # Return the complement of the regions
    # where total_region_length is the length of the whole contig the regions lie on

    # we call whatever the complement is "intron" even though depending it could be intergenic
    starts = region_array[0]
    ends = region_array[1]
    if starts[0] == 1:
        if ends[-1] == total_region_length:
            # No intron before or after the exons
            intron_starts = ends[:-1] + 1
            intron_ends = starts[1:] - 1
        else:
            # An intron after the exons but none before the exons
            intron_starts = ends + 1
            intron_ends = numpy.concatenate((starts[1:] - 1, [total_region_length]))
    else:
        if ends[-1] == total_region_length:
            # An intron before the exons but none after
            intron_starts = numpy.concatenate(([1], ends[:-1] + 1))
            intron_ends = starts - 1
        else:
            # An intron both before and after the exons
            intron_starts = numpy.concatenate(([1], ends + 1 ))
            intron_ends = numpy.concatenate((starts - 1, [total_region_length]))
    return numpy.vstack([intron_starts, intron_ends])


def annotate_regions(all_regions, all_introns, gene_info):
    # Find the introns that contain in a regionn
    # `introns` is a tuple of (ID, intron_starts, intron_ends)
    # `regions` is a region array, see merge_regions.
    # Regions must be sorted by start and pairwise disjoint
    # AND completely contained in one of the introns
    # returns a list of lists for each region listing the introns it came from
    # as (transcriptID, intron_num) pairs

    effective_lengths = collections.Counter()
    all_annotations = dict()

    for chrom in all_regions.keys():
        regions = all_regions[chrom]
        introns = all_introns[chrom]

        annotations = [[] for region in range(regions.shape[1])]
        for ID, (intron_starts, intron_ends) in introns:
            for intron_num,(start,end) in enumerate(zip(intron_starts, intron_ends)):
                # Find the first region that starts after our start
                first_after_start = numpy.searchsorted(regions[0], start, side="left") #left so we include when region start = our start,
                # Find the first region that starts after our end
                first_after_end = numpy.searchsorted(regions[0], end, side="right") # right so we exclude when region start = our end, want ones with NO overlap

                if first_after_start == first_after_end:
                    # No regions overlap our intron
                    # technically not necessary as the next step will do nothing if these are equal (the range is empty)
                    continue

                # otherwise add our transcript to annotation of every region between first and last
                [annotations[index].append((ID, intron_num)) for index in range(first_after_start, first_after_end)]

        # now go through the annotations and remove those that come from all but the most recently started gene
        # ie. if a gene is contained in another gene, the mintrons of the smaller gene should contribute only to the introns of the smaller gene
        # For mintrons from genes with no overlapping genes, this does nothing
        def remove_earlier_genes(annots):
            def get_gene_start(annot):
                (gene_id, transcript_id), intron_num = annot
                start, end = gene_info[gene_id]
                return start
            (last_gene, _), _ = max(annots, key = get_gene_start, default = ((None, None), None))
            new_annots = [((gene_id, transcript_id), intron_num) for (gene_id, transcript_id), intron_num in annots if gene_id == last_gene]
            return new_annots
        annotations = [remove_earlier_genes(annots) for annots in annotations]

        # Compute the effective lengths of each intron
        # by summing lengths of mintrons within it
        for mintron, annots in zip(regions.T, annotations):
            start, end = mintron
            mintron_length =  end - start + 1
            for ID, intron_num in annots:
                effective_lengths[ID] += mintron_length

        all_annotations[chrom] = annotations

    return all_annotations, effective_lengths

def shrink_regions(region_array, flank_size):
    # Reduce the size of each region in region_array by flank_size on each end
    # but keep at least 1 base in each region
    starts, ends = region_array
    new_starts = starts + flank_size
    new_ends = ends - flank_size

    # fix the parts where the region was less than 2*flank_size in length
    middles = numpy.floor_divide(starts + ends, 2)
    crossovers = new_starts > new_ends
    new_starts[crossovers] = middles[crossovers]
    new_ends[crossovers] = middles[crossovers]

    return numpy.vstack((new_starts, new_ends))

def flanking_regions(transcript_info, genics, flank_size):
    # Find flanking regions around our transcript that are actually around the enclosing merged genic region
    # I.e we take the region of +- flank size of our transcript
    # and then remove the full merged genic region that our transcript lives in out of that
    # and then we clip the flanking region away from other genic regions so that we don't overlap any flanks with eachother
    # and return a dictionary of {transcriptID -> ((left_start, left_end), (right_start, right_end))}
    # mapping transcript id's to left and right flanks
    flanks = dict()
    genic_starts, genic_ends = genics
    for ID, start, end in transcript_info:
        next_genic_region = numpy.searchsorted(genic_starts, start, side="right")
        our_genic_region = next_genic_region - 1
        if our_genic_region > 0:
            # Flanking region on the left
            previous_genic_region = our_genic_region - 1
            mid = numpy.floor_divide(genic_ends[previous_genic_region] + genic_starts[our_genic_region], 2)
            left_flank = (numpy.clip(mid, start - flank_size, genic_starts[our_genic_region] - 1), genic_starts[our_genic_region] - 1) # Restart..
        else: # No regions before us
            # Minimal flanking region on the left of 1 base
            left_flank = (genic_starts[our_genic_region] - flank_size - 1, genic_starts[our_genic_region] - 1)

        if next_genic_region < genics.shape[1]:
            # Flanking region on the right
            mid = numpy.floor_divide(genic_starts[next_genic_region] + genic_ends[our_genic_region], 2)
            right_flank = (genic_ends[our_genic_region] + 1, numpy.clip(genic_ends[our_genic_region] + 1, end + flank_size, mid)) # Restart..
        else: # No regions after us
            right_flank = (genic_ends[our_genic_region] + 1 , genic_ends[our_genic_region] + 1 + flank_size)

        flanks[ID] = (left_flank, right_flank)
    return flanks

if __name__ == '__main__':
    geneinfo_file_name = "/project/itmatlab/for_cris/Baby.Test_mouse_samples/resources/index_files/baby_genome.mm9/annotations"
    bamfile = "/project/itmatlab/for_cris/Baby.Test_mouse_samples/data/expression/baby_genome.mm9/aligned/baby_sample1.bam"

    intron_quant = IntronQuantificationStep(geneinfo_file_name)
    intron_quant.execute(bamfile, "", forward_read_is_sense = False)
