import collections
import itertools
import os
import argparse

import numpy
import pysam

OUTPUT_INTRON_FILE_NAME = "intron_quantifications.txt"
OUTPUT_INTRON_ANTISENSE_FILE_NAME = "intron_antisense_quantifications.txt"
OUTPUT_INTERGENIC_FILE_NAME = "intergenic_quantifications.txt"


class IntronQuantificationStep:
    def __init__(self, geneinfo_file_path, flank_size=1500):
        """
        geneinfo_file_path is a path to a tab-separated file with the following fields:
        chrom  strand  txStart  txEnd  exonCount  exonStarts  exonEnds  transcriptID  geneID  geneSymbol  biotype
        """
        self.info = None
        self.intron_normalized_antisense_counts = collections.Counter()
        self.intron_normalized_counts = collections.Counter()
        self.transcript_intron_antisense_counts = collections.Counter()
        self.transcript_intron_counts = collections.Counter()
        self.geneinfo_file_path = geneinfo_file_path
        self.flank_size = flank_size
        self.intergenic_read_counts = collections.defaultdict(collections.Counter)

    def execute(self, aligned_file_path, output_directory, forward_read_is_sense):
        # (chrom, strand) -> {(transcriptID, intron_number) -> read_count}
        intron_read_counts = collections.Counter()
        intron_antisense_read_counts = collections.Counter()
        # chrom -> {intron -> read count}

        with pysam.AlignmentFile(aligned_file_path, "rb") as alignments:

            chrom_lengths = dict(zip(alignments.references, alignments.lengths))

            #  Read in the annotation information
            self.info = AnnotationInfo(self.geneinfo_file_path, chrom_lengths)

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
                if chrom not in self.info.intergenics:
                    if chrom not in skipped_chromosomes:
                        print(f"Alignment from chromosome {chrom} skipped")
                        skipped_chromosomes.append(chrom)
                    continue

                # According to the SAM file specification, this CAN fail but I don't understand why it would
                # so just throw this assert in to verify that it doesn't, at least for now
                assert read.is_reverse != mate.is_reverse

                # Figure out the fragment's strand - depends on whether the forward or reverse reads are 'sense'
                read1_reverse_aligned = (read.is_reverse and read.is_read1) or (not read.is_reverse and read.is_read2)
                if forward_read_is_sense:
                    strand = "-" if read1_reverse_aligned else "+"
                else:
                    strand = "+" if read1_reverse_aligned else "-"
                antisense_strand = "-" if strand is "+" else "+"

                mintron_starts, mintron_ends = self.info.mintron_extents_by_chrom[chrom, strand]
                antisense_mintron_starts, antisense_mintron_ends = self.info.mintron_extents_by_chrom[chrom, antisense_strand]
                intergenic_starts, intergenic_ends = self.info.intergenic_extents_by_chrom[chrom]

                mintrons_touched = set()
                antisense_mintrons_touched = set()
                intergenics_touched = set()

                # check what regions our blocks touch
                for start, end in itertools.chain(read.get_blocks(), mate.get_blocks()):
                    # NOTE: pysam is working in 0-based, half-open coordinates! We are using 1-based
                    start = start + 1

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
                        mintrons_touched.update(range(last_mintron_before + 1, first_mintron_after))

                    # Now do the same thing as above for antisense regions
                    last_antisense_mintron_before = numpy.searchsorted(antisense_mintron_starts, start, side="right") - 1
                    first_antisense_mintron_after = numpy.searchsorted(antisense_mintron_starts, end, side="right")
                    if last_antisense_mintron_before == -1:
                        antisense_mintrons_touched.update(range(0, first_antisense_mintron_after))
                    elif antisense_mintron_ends[last_antisense_mintron_before] >= start:
                        antisense_mintrons_touched.update(range(last_antisense_mintron_before, first_antisense_mintron_after))
                    else:
                        antisense_mintrons_touched.update(range(last_antisense_mintron_before + 1, first_antisense_mintron_after))

                    # Now do the same thing as above for intergenic regions
                    last_intergenic_before = numpy.searchsorted(intergenic_starts, start, side="right") - 1
                    first_intergenic_after = numpy.searchsorted(intergenic_starts, end, side="right")
                    if last_intergenic_before == -1:
                        intergenics_touched.update(range(0, first_intergenic_after))
                    elif intergenic_ends[last_intergenic_before] >= start:
                        intergenics_touched.update(range(last_intergenic_before, first_intergenic_after))
                    else:
                        intergenics_touched.update(range(last_intergenic_before + 1, first_intergenic_after))

                # Accumulate the reads
                for mintron_index in mintrons_touched:
                    for intron in self.info.mintrons_by_chrom[chrom, strand][mintron_index].primary_introns:
                        intron_read_counts[intron] += 1

                for antisense_mintron_index in antisense_mintrons_touched:
                    for intron in self.info.mintrons_by_chrom[chrom, antisense_strand][antisense_mintron_index].primary_antisense_introns:
                        intron_antisense_read_counts[intron] += 1

                for intergenic in intergenics_touched:
                    self.intergenic_read_counts[chrom][intergenic] += 1


        # Normalize reads by effective transcript lengths
        for intron, count in intron_read_counts.items():
            effective_count = count / intron.effective_length * 1000 # FPK (fragments per kilo-base)
            self.intron_normalized_counts[intron] = effective_count
            # TODO: is this the right normalization factor for antisense reads?
            #  currently use the total length of all antisense mintrons, but this seems wrong....
            antisense_count = intron_antisense_read_counts[intron]
            if antisense_count > 0:
                self.intron_normalized_antisense_counts[intron] = antisense_count / intron.antisense_effective_length * 1000

        # Now remove counts from non-primary introns from mintrons, under the assumption that
        # if two introns overlap, we want to "subtract out" one of them from the overlap
        # We process introns in order of their transcript start: the idea being that we are modifying their
        # intron counts, so we had better be consistent as we process, going from 5' to 3' ends
        for contig, transcripts in self.info.transcripts_by_chrom.items():
            for transcript in transcripts:
                for intron in transcript.introns:
                    count = self.intron_normalized_counts.get(intron, 0)
                    if count == 0:
                        continue

                    for mintron in intron.mintrons:
                        if intron not in mintron.primary_introns:
                            for other_intron in mintron.primary_introns:
                                other_counts = self.intron_normalized_counts[other_intron]
                                # Remove the double-counts but never give negative expression
                                self.intron_normalized_counts[other_intron] = min(other_counts - count, 0)

                    # Same for anti-sense
                    antisense_count = self.intron_normalized_counts.get(intron, 0)
                    for antisense_mintron in intron.antisense_mintrons:
                        if intron not in antisense_mintron.primary_antisense_introns:
                            for other_intron in antisense_mintron.primary_antisense_introns:
                                other_counts = self.intron_normalized_antisense_counts[other_intron]
                                # Remove the double-counts but never give negative expression
                                self.intron_normalized_antisense_counts[other_intron] = min(other_counts - antisense_count, 0)

        # Transcript-level intron quantifications
        for intron, count in self.intron_normalized_counts.items():
            self.transcript_intron_counts[intron.transcript] += count
        for intron, count in self.intron_normalized_antisense_counts.items():
            self.transcript_intron_antisense_counts[intron.transcript] += count

        # Write out the results to output file
        # SENSE INTRON OUTPUT
        output_file_path = os.path.join(output_directory, OUTPUT_INTRON_FILE_NAME)
        with open(output_file_path, "w") as output_file:
            output_file.write("#gene_id\ttranscript_id\tchr\tstrand\ttranscript_intron_reads_FPK\tintron_reads_FPK\n")
            # take transcripts from all chromosomes and combine them, sorting by gene id and then transcript id
            transcript_ids = sorted((transcript.gene.gene_id, id, transcript) for id, transcript in self.info.transcripts.items())
            for (gene_id, transcript_id, transcript) in transcript_ids:

                total_count = self.transcript_intron_counts[transcript]
                intron_counts = [str(self.intron_normalized_counts[intron]) for intron in transcript.introns]


                output_file.write('\t'.join([gene_id,
                                             transcript_id,
                                             transcript.chrom,
                                             transcript.strand,
                                             str(total_count),
                                             ','.join(intron_counts),
                                             ]) + '\n')

        # ANTISENSE INTRON OUTPUT
        output_file_path = os.path.join(output_directory, OUTPUT_INTRON_ANTISENSE_FILE_NAME)
        with open(output_file_path, "w") as output_file:
            output_file.write("#gene_id\ttranscript_id\tchr\tstrand\ttranscript_intron_reads_FPK\tintron_reads_FPK\n")
            # take transcripts from all chromosomes and combine them, sorting by gene id and then transcript id
            transcript_ids = sorted((transcript.gene.gene_id, id, transcript) for id, transcript in self.info.transcripts.items())
            for (gene_id, transcript_id, transcript) in transcript_ids:

                total_count = self.transcript_intron_antisense_counts[transcript]
                intron_counts = [str(self.intron_normalized_antisense_counts[intron]) for intron in transcript.introns]


                output_file.write('\t'.join([gene_id,
                                             transcript_id,
                                             transcript.chrom,
                                             transcript.strand,
                                             str(total_count),
                                             ','.join(intron_counts),
                                             ]) + '\n')

        # TODO: do we need to normalize intergenic regions?
        #   Not clear that just dividing by their length is right since usually you have just
        #   bits and pieces expressed throughout
        output_intergenic_file_path = os.path.join(output_directory, OUTPUT_INTERGENIC_FILE_NAME)
        with open(output_intergenic_file_path, "w") as output_file:
            output_file.write("#chromosome\tintergenic_region_number\tstart\tend\treads_FPK\n")
            # take transcripts from all chromosomes and combine them, sorting by gene id and then transcript id
            chroms_sorted = sorted(self.info.intergenics.keys())
            for chrom in chroms_sorted:
                intergenics = self.info.intergenics[chrom]
                for i, intergenic in enumerate(intergenics):
                    count = self.intergenic_read_counts[chrom][i]

                    output_file.write('\t'.join([chrom,
                                                 str(i),
                                                 str(intergenic.start),
                                                 str(intergenic.end),
                                                 str(count),
                                                 ]) + '\n')


##### "Plain old data" classes representing the elements inside a gene info file

class Gene:
    """Representation of a gene, including it's transcripts.

    start, end coordinates indicate the min/max of the start/end coordinates of all its transcripts
    """
    def __init__(self, info, gene_id, chrom, strand, start, end, transcripts=None):
        if transcripts is None:
            transcripts = []
        self.info = info
        self.chrom = chrom
        self.strand = strand
        self.gene_id = gene_id
        self.start = start
        self.end = end
        self.transcripts = transcripts

    def __eq__(self, other):
        return self.gene_id == other.gene_id

    def __repr__(self):
        return f"Gene({self.gene_id}, {self.chrom}, {self.strand}, {self.start}, {self.end})"


class Transcript:
    """ Transcript of a gene, tracks its introns, exons
    """
    def __init__(self, info, gene_id, transcript_id, chrom, strand, start, end, exons=None, introns=None):
        if introns is None:
            introns = []
        if exons is None:
            exons = []
        self.info = info
        self.chrom = chrom
        self.strand = strand
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.start = start
        self.end = end
        self.exons = exons
        self.introns = introns

    def get_gene(self):
        return self.info.genes[self.gene_id]

    def __repr__(self):
        return f"Transcript({self.gene_id}, {self.transcript_id}, {self.chrom}, {self.strand}, {self.start}, {self.end})"

    gene = property(get_gene)


class Region:
    """ Any genomic span of a chromosome
    strand = +,-, or . if not strand-specific region (eg: an intergenic region)
    'comment' is any extra information to carry along, for the purposes of debugging/printing
    """
    def __init__(self, info, chrom, strand, start, end, comment=None):
        self.info = info
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.end = end
        self.comment = comment

    def __repr__(self):
        return f"{type(self).__name__}({self.chrom}, {self.strand}, {self.start}, {self.end}, {self.comment})"


class TranscriptRegion(Region):
    """ Region that is part of a transcript (eg: intron or exon)
    """
    def __init__(self, info, gene_id, transcript_id, *args):
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.effective_length = 0
        self.antisense_effective_length = 0
        Region.__init__(self, info, *args)

    def get_gene(self):
        return self.info.genes[self.gene_id]

    def get_transcript(self):
        return self.info.transcripts[self.transcript_id]

    gene = property(get_gene)
    transcript = property(get_transcript)


class Intron(TranscriptRegion):
    def __init__(self, *args):
        TranscriptRegion.__init__(self, *args)
        self.mintrons = [] # List of all mintrons that this overlaps
        self.antisense_mintrons = [] # List of mintrons on the other strand that aren't sense for something else


class Mintron(Region):
    def __init__(self, info, *args):
        Region.__init__(self, info, *args)
        self.introns = []
        self.primary_introns = []
        self.primary_gene = None
        self.primary_antisense_introns = []
        self.antisense_introns = []

    def __repr__(self):
        return f"Mintron({self.chrom}, {self.strand}, {self.start}, {self.end}, {self.primary_introns}, {self.primary_gene})"


class AnnotationInfo:
    """ Data structure containing all the information in a gene info file.

     Stores genes, transcripts, intergenic regions, and exons both for easy access by ID and for quick lookup by
     position, using sorted lists by start position, allowing binary search."""

    def __init__(self, geneinfo_file_path, chrom_lengths, flank_size=1500):
        """
        Load geneinfo file as a large datastructure containing all the information
        geneinfo file is in the bed format with 1-based, inclusive coordinates for ranges
        For example, an Ensembl gtf converted to a bed file with the provided script works.

        :param geneinfo_file_path: file path to find the geneinfo file
        :param chrom_lengths: dictionary of {chromosome -> chromosome lengths}, obtainable from header of a SAM/BAM file
        """
        self.chrom_lengths = chrom_lengths


        self.genes_by_chrom = collections.defaultdict(dict)  # All genes in a chrom
        self.genes = dict()  # All genes
        self.genic_regions_by_chrom = collections.defaultdict(collections.deque)  # genes sorted by start position
        self.transcripts = dict()  # All transcripts
        self.transcripts_by_chrom = collections.defaultdict(collections.deque)  # All transcripts in a chromosome, sorted by start position
        self.exons_by_chrom = collections.defaultdict(collections.deque)  # Exons sorted by start position
        with open(geneinfo_file_path) as geneinfo_file:
            for line in geneinfo_file:
                if line.startswith("#"):
                    continue
                chrom, strand, tx_start, tx_end, exon_count, exon_starts, exon_ends, transcript_id, gene_id, *other = line.split('\t')
                tx_start, tx_end = int(tx_start), int(tx_end)
                exon_starts = [int(start) for start in exon_starts.split(',')]
                exon_ends = [int(end) for end in exon_ends.split(',')]
                exons = [TranscriptRegion(self, gene_id, transcript_id, chrom, strand, start, end, "exon") for
                         start, end in zip(exon_starts, exon_ends)]

                intron_starts = [end + 1 for end in exon_ends[:-1]]
                intron_ends = [start - 1 for start in exon_starts[1:]]
                introns = [Intron(self, gene_id, transcript_id, chrom, strand, start, end, "intron") for
                           start, end in zip(intron_starts, intron_ends)]

                transcript = Transcript(self, gene_id, transcript_id, chrom, strand, tx_start, tx_end, exons, introns)
                self.transcripts[transcript_id] = transcript
                self.transcripts_by_chrom[chrom, strand].append(transcript)

                if gene_id in self.genes:
                    transcript.gene.start = min(transcript.gene.start, tx_start)
                    transcript.gene.end = max(transcript.gene.end, tx_end)
                    transcript.gene.transcripts.append(transcript)
                else:
                    gene = Gene(self, gene_id, chrom, strand, tx_start, tx_end, transcripts=[transcript])
                    self.genes[gene_id] = gene
                    self.genes_by_chrom[chrom, strand][gene_id] = gene
                    self.genic_regions_by_chrom[chrom, strand].append(gene)

                self.exons_by_chrom[chrom, strand].extend(exons)

        # Sort regions by start position
        self.exons_by_chrom = {contig: sorted(exons, key=lambda x: x.start)
                               for contig, exons in self.exons_by_chrom.items()}
        self.genic_regions_by_chrom = {contig: sorted(genics, key=lambda x: x.start)
                                       for contig, genics in self.genic_regions_by_chrom.items()}

        self.transcripts_by_chrom = {contig: sorted(transcripts, key=lambda x: x.start)
                                     for contig, transcripts in self.transcripts_by_chrom.items()}

        # Merge exons
        self.merged_exons = self.merge_regions(self.exons_by_chrom)
        self.merged_genics = self.merge_regions(self.genic_regions_by_chrom, merge_strands=True)

        # Compute the flanks of each gene (if any)
        # and add them to self.merged_genics and to the each Gene as a new intron (if any)
        self.flanked_genics = self.add_flanks(flank_size)

        # Get intergenic regions sorted by start
        self.intergenics = self.complement_regions(self.flanked_genics, by_strand=False)

        # Compute the mintrons and their inclusions
        # mintrons are each region that is neither intergenic nor exonic - pieces of intron
        intergenic_and_exonic = {(chrom, strand): sorted(itertools.chain(self.intergenics[chrom], self.merged_exons[chrom, strand]), key=lambda x: x.start)
                                    for chrom, strand in self.merged_exons}
        self.intergenic_and_exonic = self.merge_regions(intergenic_and_exonic)
        self.mintrons_by_chrom = self.complement_regions(self.intergenic_and_exonic, cls=Mintron)

        # The start,ends of the mintrons/intergenics as numpy arrays for fast lookup
        self.mintron_extents_by_chrom = {(chrom, strand): (numpy.array([mintron.start for mintron in mintrons]),
                                                          numpy.array([mintron.end for mintron in mintrons]))
                                            for (chrom, strand), mintrons in self.mintrons_by_chrom.items()}
        self.intergenic_extents_by_chrom = {chrom: (numpy.array([intergenic.start for intergenic in intergenics]),
                                                              numpy.array([intergenic.end for intergenic in intergenics]))
                                                for chrom, intergenics in self.intergenics.items()}

        # Give mintrons their annotations
        for (chrom, strand), transcripts in self.transcripts_by_chrom.items():
            mintrons = self.mintrons_by_chrom[chrom, strand]
            mintron_starts = numpy.array([m.start for m in mintrons])
            for transcript in transcripts:
                for intron in transcript.introns:
                    # Find corresponding mintrons
                    # first find the one that starts before us, which we may or may not intersect
                    index_before = numpy.searchsorted(mintron_starts, intron.start, side="right") - 1
                    if index_before == -1:  # No mintrons before our start
                        index_before = 0
                    mintron_before = mintrons[index_before]
                    if mintron_before.end < intron.start:
                        # We do not intersect with this one
                        # But the next one either does intersect or is past us to the right
                        index_before = index_before + 1

                    # Find first mintron that is completely to our right
                    index_after = numpy.searchsorted(mintron_starts, intron.end, side="right")
                    if index_after == len(mintron_starts):  # no mintrons after our end
                        index_after = len(mintron_starts)

                    mintrons_intersected = range(index_before,
                                                 index_after)  # possibly nothing in this range if these are equal

                    intron.mintrons = [mintrons[i] for i in mintrons_intersected]
                    for mintron in intron.mintrons:
                        mintron.introns.append(intron)
                        if mintron.primary_gene is None:
                            mintron.primary_gene = intron.gene
                            mintron.primary_introns = [intron]
                        elif mintron.primary_gene == intron.gene:
                            mintron.primary_introns.append(intron)
                        else:
                            if intron.gene.start > mintron.primary_gene.start:
                                # Overwrite existing primary introns since we start closer to their mintron
                                mintron.primary_introns = [intron]
                                mintron.primary_gene = intron.gene
                            else:
                                pass
                                # Skip if our gene starts before the primary one
                                # Idea is if a gene lives inside another one, it should get priority over the introns
                                # of the bigger gene since we won't want to count a read for both genes

            # Now add in anti-sense if necessary
            # Sense gets preference so annotate only if no sense annotation already exists
            antisense_strand = "+" if strand == "-" else "-"
            for transcript in self.transcripts_by_chrom[chrom, antisense_strand]:
                for intron in transcript.introns:
                    # Find corresponding mintrons
                    # first find the one that starts before us, which we may or may not intersect
                    index_before = numpy.searchsorted(mintron_starts, intron.start, side="right") - 1
                    if index_before == -1:  # No mintrons before our start
                        index_before = 0
                    mintron_before = mintrons[index_before]
                    if mintron_before.end < intron.start:
                        # We do not intersect with this one
                        # But the next one either does intersect or is past us to the right
                        index_before = index_before + 1

                    # Find first mintron that is completely to our right
                    index_after = numpy.searchsorted(mintron_starts, intron.end, side="right")
                    if index_after == len(mintron_starts):  # no mintrons after our end
                        index_after = len(mintron_starts)

                    mintrons_intersected = range(index_before,
                                                 index_after)  # possibly nothing in this range if these are equal

                    intron.antisense_mintrons = [mintrons[i] for i in mintrons_intersected]
                    for mintron in intron.antisense_mintrons:
                        mintron.antisense_introns.append(intron)
                        if mintron.primary_gene is None:
                            mintron.primary_antisense_introns = [intron]
                            mintron.primary_gene = intron.gene
                        elif mintron.primary_gene == intron.gene:
                            mintron.primary_antisense_introns.append(intron)
                        elif mintron.primary_gene.strand == antisense_strand:
                            if intron.gene.start > mintron.primary_gene.start:
                                # Overwrite existing primary introns since we start closer to their mintron
                                mintron.primary_antisense_introns = [intron]
                                mintron.primary_gene = intron.gene
                            else:
                                pass
                                # Skip if our gene starts before the primary one
                                # Idea is if a gene lives inside another one, it should get priority over the introns
                                # of the bigger gene since we won't want to count a read for both genes
                        else:
                            pass
                            # Do nothing since the mintron already has sense annotations, don't give it any anti-sense

            # Compute effective lengths of each intron
            for transcript in self.transcripts.values():
                for intron in transcript.introns:
                    intron.effective_length = 0
                    for mintron in intron.mintrons:
                        intron.effective_length += mintron.end - mintron.start + 1
                    intron.antisense_effective_length = 0
                    for mintron in intron.antisense_mintrons:
                        intron.antisense_effective_length += mintron.end - mintron.start + 1



    def add_flanks(self, flank_size):
        """
        Add flanks to each genic region up to size flank_size on each end
        These account for reads that go past the "ends" of the gene but should be thought of as belonging to that gene
        rather than to intergenic regions.
        Flanks are added as first/last introns to every transcript, so all transcripts will have at least two introns.
        Sometimes there is no room for a flank to be added, in which case the flank is given start/stop coordinates
        of 0,0.

        :param flank_size: Size of flanks to add
        :return: flanked_genics, dictionary of genic regions, sorted by start position with flanks added
        """
        # Look through each gene, check if contained in other genes on either side
        # on sides not contained, add a flank of size flank_size
        # up to half the size of the gap until the next exon

        flanked_genics_by_chrom = dict()
        for (chrom, strand), genes in self.genes_by_chrom.items():
            genic_regions = self.genic_regions_by_chrom[chrom, strand]
            starts = numpy.array([g.start for g in genic_regions])
            for gene in genes.values():
                # Find genic region that we are included in
                index = numpy.searchsorted(starts, gene.start,
                                           side="right") - 1  # Always >= 0 since a gene's own genic region starts before the gene
                our_genic_region = genic_regions[index]

                # Left flank
                cmt = "left flank"
                if our_genic_region.start != gene.start:
                    # no flank on left, since some other gene covers the left
                    flank_start = 0
                    flank_end = 0
                else:
                    # Genic region on the left, extend to half that distance
                    if index > 0:
                        previous_genic_region = genic_regions[index - 1]
                        mid = (previous_genic_region.end + our_genic_region.start) // 2 + 1 # midpoint (rounded down)
                        flank_start = max(mid, gene.start - flank_size)
                    else:
                        flank_start = max(1, gene.start - flank_size)
                    flank_end = our_genic_region.start - 1
                    if flank_start == flank_end:
                        flank_start, flank_end = 0, 0  # no flank
                        cmt = "no left flank"
                for transcript in gene.transcripts:
                    left_flank = Intron(self, transcript.gene_id, transcript.transcript_id, chrom, strand,
                                                  flank_start, flank_end, cmt)
                    transcript.introns.insert(0, left_flank)

                # Right flank
                cmt = "right flank"
                if our_genic_region.end != gene.end:
                    # no flank on right, since some other gene covers the right
                    flank_start = 0
                    flank_end = 0
                else:
                    # Genic region on the right, extend to half that distance
                    if index < len(genic_regions) - 1:
                        next_genic_region = genic_regions[index + 1]
                        # Always leave a gap of at least 1 bases by an extra + 1 here
                        mid = (next_genic_region.start + our_genic_region.end) // 2 - 1 # midpoint (rounded up)
                        flank_end = min(mid, gene.end + flank_size)
                    else:
                        flank_end = min(self.chrom_lengths[chrom], gene.end + flank_size)
                    flank_start = our_genic_region.end + 1
                    if flank_start == flank_end:
                        flank_start, flank_end = 0, 0  # no flank
                        cmt = "no right flank"
                for transcript in gene.transcripts:
                    right_flank = Intron(self, transcript.gene_id, transcript.transcript_id, chrom, strand,
                                                   flank_start, flank_end, cmt)
                    transcript.introns.append(right_flank)

            # Now repeat some of that to update the merged genic regions with their flanks
            flanked_genics = collections.deque()
            for index, region in enumerate(genic_regions):
                if index == 0:  # no previous genic region
                    left_start = max(1, region.start - flank_size)
                else:
                    previous_region = genic_regions[index - 1]
                    mid = (previous_region.end + region.start) // 2 + 1
                    left_start = max(mid, region.start - flank_size)

                if index == len(genic_regions) - 1:  # no next genic region
                    right_end = min(self.chrom_lengths[chrom], region.end + flank_size)
                else:
                    next_region = genic_regions[index + 1]
                    mid = (region.end + next_region.start) // 2 - 1
                    right_end = min(mid, region.end + flank_size)
                flanked_genics.append(Region(self, chrom, ".", left_start, right_end, comment="flanked genic"))

            flanked_genics_by_chrom[chrom] = flanked_genics
        return flanked_genics_by_chrom

    def merge_regions(self, regions_by_chrom, merge_strands=False):
        merged_by_chrom = dict()
        for (chrom, strand), regions in regions_by_chrom.items():
            # merge regions in this contig
            last_start = -1
            last_end = -1
            merged = collections.deque()
            merged_regions = None
            for region in regions:
                start, end = region.start, region.end
                if start <= last_end + 1:
                    # Overlaps with previous start...
                    # ... OR is adjacent to previous exon (want to join touching exons)
                    if end > last_end:
                        # but extends further
                        merged_regions.append(region)
                        merged[-1] = Region(self, chrom, strand, last_start, end, comment=merged_regions)
                        last_end = end
                    else:
                        pass  # Skip, contained in the existing exon
                else:
                    # New exon, no overlap
                    merged_regions = [region]
                    merged.append(Region(self, chrom, strand, start, end, comment=merged_regions))
                    last_start = start
                    last_end = end

            if merge_strands:
                merged_by_chrom[chrom] = merged
            else:
                merged_by_chrom[chrom, strand] = merged
        return merged_by_chrom

    def complement_regions(self, regions_by_chrom, by_strand=True, cls=Region):
        # Assumes regions_by_chrom are merged
        complements_by_chrom = dict()
        for contig, regions in regions_by_chrom.items():
            if by_strand:
                chrom, strand = contig
            else:
                chrom = contig
                strand = "."

            complements = collections.deque()
            start = 1
            for region in regions:
                if region.start == 1:
                    # No complement on left end
                    start = region.end + 1
                    continue
                end = region.start - 1
                # Never should get start == end here since we are assuming that adjecent regions have been merged
                complements.append(cls(self, chrom, strand, start, end))
                start = region.end + 1

            if start < self.chrom_lengths[chrom]:
                # need one more region on the end
                complements.append(cls(self, chrom, strand, start, self.chrom_lengths[chrom]))

            if by_strand:
                complements_by_chrom[chrom, strand] = complements
            else:
                complements_by_chrom[chrom] = complements
        return complements_by_chrom


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam_file", help="BAM or SAM file of strand-specific genomic alignment")
    parser.add_argument("-i", "--info_file", help="Geneinfo file in BED format with 1-based, inclusive coordinates")
    parser.add_argument("-o", "--output_directory", help=f"Directory where to output files {OUTPUT_INTRON_FILE_NAME}, {OUTPUT_INTRON_ANTISENSE_FILE_NAME}, {OUTPUT_INTERGENIC_FILE_NAME}")
    parser.add_argument("--forward_read_is_sense", help="Set if forward read is sense. Default is False. Strand-specificity is assumed", action="store_const", const=True, default=False)

    args = parser.parse_args()

    intron_quant = IntronQuantificationStep(args.info_file)
    intron_quant.execute(args.bam_file, args.output_directory, args.forward_read_is_sense)

"""
Example:
python intron_quant.py -b "/Users/tgbrooks/BEERS2.0/data/baby_sample1.bam" -i "/Users/tgbrooks/BEERS2.0/data/baby_genome.mm9/annotations" -o .
"""
