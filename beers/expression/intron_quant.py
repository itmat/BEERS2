import collections
import pysam
import numpy

class IntronQuantificationStep:
    def __init__(self, geneinfo_file_path):
        '''
        geneinfo_file_path is a path to a tab-separated file with the following fields:
        chrom  strand  txStart  txEnd  exonCount  exonStarts  exonEnds  transcriptID  geneID  geneSymbol  biotype
        '''
        self.geneinfo_file_path = geneinfo_file_path

    def execute(self, aligned_file_path, output_directory):
        nonexonic_counts = {chrom: collections.defaultdict(lambda: 0) for chrom in self.exons}

        with pysam.AlignmentFile(aligned_file_path, "rb") as alignments:

            # Read in the gene alignment information
            chrom_lengths = dict(zip(alignments.references, alignments.lengths))
            exons, mintrons, genics, intergenics, mintron_annotations = load_gene_info(self.geneinfo_file_path, chrom_lengths)

            # (chrom, strand) -> {(transcriptID, intron_number) -> read_count}
            intron_read_counts = collections.defaultdict(lambda: collections.defaultdict(lambda:0))
            # chrom -> {intron -> read count}
            intergenic_read_counts = collections.defaultdict(lambda: collections.defaultdict(lambda:0))

            # Go through all reads, and compare
            skipped_chromosomes = []
            for read in alignments.fetch(until_eof=True):
                # Use mapped, uniquely mapped reads only
                if read.is_unmapped() and read.get_tag("NH") == 1:
                    continue

                chrom = read.reference_name()
                if chrom not in self.nonexons: # Check annotations are available for that chromosome
                    if chrom not in skipped_chromosomes:
                        print(f"Alignment from chromosome {chrom} skipped")
                        skipped_chromosomes.append(chrom)
                    continue

                strand = "-" if read.is_reverse() else "+"

                mintron_starts, mintron_ends = mintrons[(chrom, strand)]

                mintrons_touched = set()
                intergenics_touched = set()

                # check what regions our blocks touch
                for start,end in read.get_blocks:
                    # NOTE: pysam is working in 0-based, half-open coordinates! We are using 1-based
                    start = start+1

                    # We may intersect this mintron, depending on where its end lies
                    last_mintron_before = numpy.searchsorted(mintron_starts, start, side="right")

                    # We definitely do not intersect this mintron, it starts after our end
                    first_mintron_after = numpy.searchsorted(mintron_starts, end, side="right") + 1
                    # but all between the two, we do intersect

                    if mintron_ends[last_mintron_before] >= start:
                        # We do intersect the last one to start before us
                        mintrons_touched.update(range(last_mintron_before, first_mintron_after))
                    else:
                        # We only start intersecting the first one after that
                        mintrons_touched.update(range(last_mintron_before+1, first_mintron_after))


                    # Now do the same thing for intergenic regions
                    last_intergenic_before = numpy.searchsorted(intergenic_starts, start, side="right")
                    first_intergenic_after = numpy.searchsorted(intergenic_starts, end, side="right") + 1
                    if intergenic_ends[last_intergenic_before] >= start:
                        intergenics_touched.update(range(last_intergenic_before, first_intergenic_after))
                    else:
                        intergenics_touched.update(range(last_intergenic_before+1, first_intergenic_after))

                # Accumulate the reads
                for mintron in mintrons_touched:
                    for intron in mintron_annotations[mintron]:
                        intron_read_counts[intron] += 1

                for intergenic in intergenics_touch:
                    intergenic_read_counts[intergenic] += 1


        return intron_read_counts, intergenic_read_counts

def load_gene_info(geneinfo_file_path, chrom_lengths):
    # Assumes that geneinfo_file is sorted by feature starts within a chromosome!!
    # and that the chromosomes in geneinfo_file are all present in chrom_lengths (i.e. in the bam file header)

    # Gather exon boundaries from all annotations, sorted by their start
    # Collect exons by chromosome and strand, and intergenics just by chromosome
    make_deque = lambda: collections.deque()
    exons = collections.defaultdict(make_deque)
    genic_regions = collections.defaultdict(make_deque)
    all_introns = collections.defaultdict(make_deque)
    with open(geneinfo_file_path) as geneinfo_file:
        for line in geneinfo_file:
            if line.startswith("#"):
                continue
            chrom, strand, txStart, txEnd, exonCount, exonStarts, exonEnds, transcriptID, geneID, geneSymbol, biotype = line.split('\t')
            txStart, txEnd = int(txStart), int(txEnd)
            exonStarts = [int(start) for start in exonStarts]
            exonEnds = [int(end) for end in exonEnds]
            exons[(chrom,strand)].extend(zip(exonEnds, exonStarts)) # Track strandedness of exons
            genic_regions.append( (txStart, txEnd) )

            intronStarts = [end + 1 for end in exonEnds[1:]]
            intronEnds = [start - 1 for start in  exonStarts[:-1]]
            all_introns[(chrom,strand)].append( (transcriptID, intronStarts, intronEnds) )

    # Convert to arrays (array[0] is array of starts, array[1] array of ends
    exons = {chrom: numpy.sort(numpy.array(exons).T) for chrom,exons in merged_exons.items()}
    genics = {chrom: numpy.sort(numpy.array(genic).T) for chrom,genics in merged_genics.items()}

    # Merge regions that overlap
    merged_exons = {chrom: merge_regions(exons) for chrom, exons in all_exons.items()}
    merged_genic = {chrom: merge_regions(genic) for chrom, genic in genic_exons.items()}

    # Intergenics are complements of all the genics
    intergenics = {chrom: complement_regions(genic, chrom_lengths[chrom]) for chrom,genics in merged_genics.items()}

    # intronic regions are ones that between these merged exons that ARENT intergenic
    # so combine exons and intergenics
    # These regions are disjoint, so just need to sort them together
    intergenics_or_exons = {(chrom,strand): numpy.concatenate(exons, intergenics[chrom]) for (chrom,strand),exons in merged_exons.items()}
    [intergenic_or_exon.sort() for intergenic_or_exon  in intergenic_or_exons.values()] # And sort them into the right order

    # We define mintrons to be a region that is not intergenic and is not in an exon of THIS strand
    # i.e. is in at least one gene on this strand and not in any exon of any gene of this strand
    # (for mini-introns, of which there might be several in a real intron since exons of other genes could break up the intron)
    # Analogously, we should define the 'exons' above to be 'maxons' since they're max of all touching exons...
    mintrons = {chrom: complement_regions(regions, chrom_lengths[chrom]) for (chrom,strand),regions in intergenic_or_exons.items()}

    # Annotate the mintrons by the transcriptID + introns
    mintron_annotations = {chrom: annotate_regions(mintrons[chrom], all_introns[chrom]) for chrom in mintrons}


    return merged_exons, mintrons, merged_genics, intergenics, mintron_annotations

def merge_regions(region_arrays):
    # Given region array, i.e. an numpy array such that
    # region_array[0] lists the start positions
    # region_array[1] lists the end positions
    # Sorted by start poisition
    # Merge any regions that overlap (or are adjacent)
    last_start = 0
    last_end = 0
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
                # Skip, contained in the existing exon
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
    starts = array[0]
    ends = array[1]
    if exon_starts[0] == 1:
        if exon_ends[-1] == chrom_length:
            # No nonexons before or after the exons
            intron_starts = exon_ends[:-1] + 1
            intron_ends = exon_starts[1:] - 1
        else:
            # A nonexon  after the exons but none before the exons
            intron_starts = exon_ends + 1
            intron_ends = numpy.concatenate(exon_starts[1:] - 1, [chrom_length])
    else:
        if exon_ends[-1] == chrom_length:
            # A nonexon  before the exons but none after
            intron_starts = exon_ends[:-1] + 1
            intron_ends = exon_starts[1:] - 1
        else:
            # A nonexon both before and after the exons (typical case, really the only case in real genomes)
            intron_starts = numpy.concatenate([0], exon_ends) + 1
            intron_ends = numpy.concatenate(exon_starts - 1, [chrom_length])
    return numpy.hstack([intron_starts, intron_ends])


def annotate_regions(regions, introns):
    # Find the introns that contain in a regionn
    # `introns` is a tuple of (ID, intron_starts, intron_ends)
    # `regions` is a region array, see merge_regions.
    # Regions must be sorted by start and pairwise disjoint
    # AND completely contained in one of the introns
    # returns a list of lists for each region listing the introns it came from
    # as (transcriptID, intron_num) pairs
    annotations = [[] for region in regions]
    for ID, intron_starts, intron_ends in introns:
        for intron_num,(start,end) in enumerate(zip(intron_starts, intron_ends)):
            # Find the first region that starts after our start
            first_after_start = numpy.searchsorted(regions, start, side="right") + 1
            # Find the first region that starts after our end
            first_after_end = numpy.searchsorted(regions, end, side="right") + 1

            if first_after_start == first_after_end:
                # No regions overlap our intron
                # technically not necessary as the next step will do nothing if these are equal (the range is empty)
                continue

            # otherwise add our transcript to annotation of every region between first and last
            [annotations[index].append((ID, intron_num)) for index in range(first_after_start, first_after_end)]
    return annotations
