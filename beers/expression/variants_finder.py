import sys
import argparse
import os
import re
from timeit import default_timer as timer
from collections import namedtuple
from operator import attrgetter, itemgetter
import math
from io import StringIO


Read = namedtuple('Read', ['type', 'chromosome', 'position', 'description'])
"""
A named tuple that possesses all the attributes of a variant
type:  match (M), deletion (D), insertion (I)
chromosome: chrN
position: position on ref genome
description: description of the variant (e.g., C, IAA, D5, etc.) 
"""


class VariantsFinder:
    """
    This class creates a text file listing variants for those locations in the reference genome having variants.
    The variants include snps and indels with the number of reads attribute to each variant.
    The text-based input file has no header and the following columns:
    1) CHROMOSOME (column 1 in a SAM file)
    2) START (column 3 in a SAM file)
    3) CIGAR  (column 4 in a SAM file)
    4) SEQ  (column 10 in a SAM file)
    The reads must be sorted by location

    This script outputs a file that gives the full breakdown at each
    location in the genome of the number of A's, C's, G's and T's as
    well as the number of each size of insertion and deletion.
    If it's an insertion the sequence of the insertion itself is given.
    So for example a line of output like the following means
    29 reads had a C in that location, one had a T and
    also one read had an insertion TT and three reads had an insertion TTT
    chr1:10128503 | C:29 | T:1 | IT:1 ITTT:3
    """

    def __init__(self, alignment_map_filename, reference_genome_filename, sort_by_entropy, depth_cutoff):
        self.entropy_sort = True if sort_by_entropy else False
        self.depth_cutoff = depth_cutoff or 10
        self.alignment_map_filename = alignment_map_filename
        self.reference_genome_filename = reference_genome_filename
        self.variants_filename = os.path.splitext(alignment_map_filename)[0] + "_variants.txt"
        self.clip_at_start_pattern = re.compile("(^\d+)[SH]")
        self.clip_at_end_pattern = re.compile("\d+[SH]$")
        self.variant_pattern = re.compile("(\d+)([NMID])")
        self.indel_pattern = re.compile("\|([^|]+)")
        self.reference_sequence_generator = self.find_reference_sequence()
        try:
            os.remove(self.variants_filename)
        except OSError:
            pass

    def remove_clips(self, cigar, seq):
        """
        Remove soft and hard clips at the beginning and end of the cigar string and remove soft and hard clips at
        the beginning of the seq as well.  Modified cigar string and sequence are returned
        :param cigar: raw cigar string from read
        :param seq: raw sequence string from read
        :return: tuple of modified cigar and sequence strings (sans clips)
        """
        clip_at_start = re.search(self.clip_at_start_pattern, cigar)
        if clip_at_start:
            cigar = re.sub(self.clip_at_start_pattern, "", cigar)
            seq = seq[int(clip_at_start.group(1)):]
        cigar = re.sub(self.clip_at_end_pattern, "", cigar)
        return cigar, seq

    def dump_to_file(self, chromosome, reads):
        """
        Parses the reads dictionary (read named tuple:read count) for each chromosome - position to create
        a line with the variants and their counts delimited by pipes.  Dumping each chromosome's worth of
        data at a time is done to avoid too sizable a dictionary.  Additionally, if the user requests a sort by entropy,
        this function will do that ordering and send that data to stdout.
        :param chromosome: chromosome whose variants are being dumped to file
        :param reads: dictionary of reads to read counts
        """

        print(f"Dumping chromosome {chromosome} reads.")

        # Keep retrieving reference data from the reference genome file until the seq for the chromosome matching the
        # one under study here is found.
        while True:
            reference_chromosome, reference_sequence = self.reference_sequence_generator.__next__()
            if reference_chromosome == chromosome:
                break

        # Open the file that will hold the variants discovered.
        with open(self.variants_filename, "a") as variants_file:

            # Initializing the variable that will hold position information objects
            position_info = None

            # This dictionary is only used if the user requests that the read lines be sorted by entropy
            entropy_map = dict()

            # Iterate over the reads in the dictionary of variants to read counts sorted by the read position
            for read in sorted(reads.keys(), key=attrgetter('position')):

                # Initial iteration - set up position information object.
                if not position_info:
                    position_info = PositionInfo(read.chromosome, read.position)

                # If the new position differs from the position of the position information currently being
                # consolidated, dump the current position information to the variants file if it is determined to
                # contain at least one variant.  In either case, create a new position information object for the new
                # position.
                if read.position != position_info.position:

                    reference_base = reference_sequence[position_info.position - 1]
                    if position_info.has_variant(reference_base):
                        variants_file.write(position_info.__str__())

                    # If the sort by entropy option is selected, also add to the entropy map dictionary the position
                    # information entropy, keyed by the line content but only if the total number of reads exceeds the
                    # depth cutoff.
                    if self.entropy_sort and position_info.get_total_reads() >= int(self.depth_cutoff):
                            entropy_map[position_info.__str__()] = position_info.calculate_entropy()

                    position_info = PositionInfo(read.chromosome, read.position)

                # Add the read description and read count to the position information
                position_info.add_read(read.description, reads[read])

        # Now that the reads are exhausted for this chromosome, dump the data from the current position information
        # object to the file.
        reference_base = reference_sequence[position_info.position - 1]
        if position_info.has_variant(reference_base):
            variants_file.write(position_info.__str__())

        # If the user selected the sort by entropy option, other the entropy_map entries in descending order
        # of entropy and print to std out.
        if self.entropy_sort:
            sorted_entropies = sorted(entropy_map.items(), key=itemgetter(1), reverse=True)
            for key, value in sorted_entropies:
                print(key, end='')

    def find_reference_sequence(self):
        """
        Returns a generator object that gets the next chromosome, sequence pair as a tuple from the reference genome.
        If the file is exhausted, the process stops with a warning suggesting that the alignment file either has reads
        for chromosomes not listed in the reference genome file or that the chromosome ordering of the two files do
        not match.
        :return: generator object provided next reference chromosome, sequence pair
        """

        with open(self.reference_genome_filename, 'r') as reference_genome_file:
            while True:
                line = reference_genome_file.readline()
                if not line:
                    break
                if line.startswith(">"):
                    chromosome = line[1:].rstrip("\n")
                    yield chromosome, reference_genome_file.readline()
        print("Warning:  Variants file not exhausted but reference genome file is.  Could be due to mis-ordering.")
        sys.exit(1)

    def find_variants(self):
        """
        Iterate over the input txt file containing cigar, seq, start location, chromosome for each read and consolidate
        reads for each position on the genome.
        """

        reads, current_chromosome = dict(), ''
        with open(self.alignment_map_filename, "r") as alignment_map_file:
            for line in alignment_map_file:
                (chromosome, start, cigar, seq) = line.rstrip('\n').split('\t')

                # Initial iteration - set up the current chromosome
                if not current_chromosome:
                    current_chromosome = chromosome

                # If we are starting a new chromosome, dump the reads for the existing chromosome into the output
                # file and renew the reads dictionary for the new chromosome.
                if current_chromosome and current_chromosome != chromosome:
                    self.dump_to_file(current_chromosome, reads)
                    reads.clear()
                    current_chromosome = chromosome

                cigar, seq = self.remove_clips(cigar, seq)
                current_pos_in_genome = int(start)
                loc_on_read = 1

                # Iterate over the variant types and lengths in the cigar string
                for match in re.finditer(self.variant_pattern, cigar):
                    length = int(match.group(1))
                    read_type = match.group(2)

                    # Skip over N type reads since these generally represent a read bracketing an intron
                    if read_type == "N":
                        current_pos_in_genome += length
                        continue

                    # For a match, record all the snps at the each location continuously covered by this read type
                    if read_type == "M":
                        stop = current_pos_in_genome + length
                        while current_pos_in_genome < stop:
                            location = current_pos_in_genome
                            reads[Read(read_type, current_chromosome, location, seq[loc_on_read - 1])] = \
                                reads.get(
                                    Read(read_type, current_chromosome, location, seq[loc_on_read - 1]), 0) + 1
                            loc_on_read += 1
                            current_pos_in_genome += 1
                        continue

                    # For a deletion, designate the read named tuple description with a Dn where n is the
                    # length of the deletion starting at this position.  In this way, subsequent reads having a
                    # deletion of the same length at the same position will be added to this key.
                    if read_type == "D":
                        location = current_pos_in_genome
                        reads[Read(read_type, current_chromosome, location, f'D{length}')] = \
                            reads.get(Read(read_type, current_chromosome, location, f'D{length}'), 0) + 1
                        current_pos_in_genome += length
                        continue

                    # For an insert, designate the read named tuple description with an Ib+ where b+ are the
                    # bases to a inserted starting with this position.  In this way, subsequent reads having an
                    # insertion of the same bases at the same position will be added to this key.
                    if read_type == "I":
                        location = current_pos_in_genome
                        insertion_sequence = seq[loc_on_read - 1: loc_on_read - 1 + length]
                        reads[Read(read_type, current_chromosome, location, f'I{insertion_sequence}')] = \
                            reads.get(
                                Read(read_type, current_chromosome, location, f'I{insertion_sequence}'), 0) + 1
                        loc_on_read += length

            # For each chromosome completed dump the dictionary contents to avoid a large memory footprint.
            self.dump_to_file(current_chromosome, reads)

    @staticmethod
    def main():
        """
        CLI Entry point into the variants_finder program.  Parses the use input, created the VariantsFinder object,
        passing in the arguments and runs the process to find the variants inside a timer.
        """
        parser = argparse.ArgumentParser(description='Make Variants File')
        parser.add_argument('-a', '--alignment_map_filename',
                            help="Textfile providing chromosome, start postion, cigar,"
                                 " and sequence only for each read.")
        parser.add_argument('-g', '--reference_genome_filename',
                            help="Path to the related reference genome.  Used to eliminate read positions that contain"
                                 "no variants.")
        parser.add_argument('-s', '--sort_by_entropy', action='store_true',
                            help="Optional request to sort line in order of descreasing entropy.")
        parser.add_argument('-c', '--cutoff_depth', type=int, default=10,
                            help="Integer to indicate minimum read depth a position must have for inclusion."
                                 " If the option is not selected, a default of 10 will be applied as the minimum"
                                 " read depth.  Note that this option is used only if the sort_by_entropy option is"
                                 " invoked.")
        args = parser.parse_args()
        print(args)
        variants_finder = VariantsFinder(args.alignment_map_filename,
                                         args.reference_genome_filename,
                                         args.sort_by_entropy,
                                         args.cutoff_depth)
        start = timer()
        variants_finder.find_variants()
        end = timer()
        sys.stderr.write(f"Variants Finder: {end - start} sec\n")


class PositionInfo:
    """
    This class is meant to capture all the read data associated with a particular chromsome and position on the
    genome.  It is used to ascertain whether this position actually holds a variant.  If it does, the data is
    formatted into a string to be written into the variants file.
    """

    def __init__(self, chromosome, position):
        self.chromosome = chromosome
        self.position = position
        self.reads = []

    def add_read(self, description, read_count):
        self.reads.append((description, read_count))

    def get_total_reads(self):
        return sum([read[1] for read in self.reads])

    def get_abundances(self):
        return [read[1] / self.get_total_reads() for read in self.reads]

    def calculate_entropy(self):
        """
        Use the top two abundances (if two) of the variants for the given position to compute an entropy.  If
        only one abundance is given, return 0.
        :return: entropy for the given position
        """
        abundances = self.get_abundances()
        if len(abundances) < 2:
            return 0

        # cloning the abundances list since lists are mutable.
        abundances_copy = abundances.copy()

        # Retrieve the highest abundance, then remove it and retrieve the highest abundance again to get the
        # next highest.
        max_abundances = [max(abundances_copy)]
        abundances_copy.remove(max_abundances[0])
        max_abundances.append(max(abundances_copy))

        # In the event that the second abundance reads nearly 0, just return 0
        if max_abundances[1] < 0.0001:
            return 0

        # Use a scale factor to normalize to the two abundances used to calculate entropy
        scale = 1 / sum(max_abundances)
        max_abundances = [scale * max_abundance for max_abundance in max_abundances]
        return -1 * max_abundances[0] * math.log2(max_abundances[0]) - max_abundances[1] * math.log2(max_abundances[1])

    def has_variant(self, reference_base):
        """
        To have a variant, the position information must contain a single read description and that description may
        not be identical to the base at that position in the reference genome.
        :param reference_base: base at this position in the reference genome.
        :return: True if the position information included at least one variant and false otherwise.
        """
        if len(self.reads) > 1 or self.reads[0][0] != reference_base:
            return True
        return False

    def __str__(self):
        """
        Provides a string representation that may be used to dump to a file.
        :return: string representation
        """
        abundances = [str(round(abundance, 2)) for abundance in self.get_abundances()]
        s = StringIO()
        s.write(f'{self.chromosome}:{self.position}')
        for read in self.reads:
            s.write(f' | {read[0]}:{read[1]}')
        s.write(f"\tTOT={self.get_total_reads()}")
        s.write(f"\t{','.join(abundances)}")
        s.write(f"\tE={self.calculate_entropy()}\n")
        return s.getvalue()


if __name__ == "__main__":
    sys.exit(VariantsFinder.main())

'''Example call
python variants_finder.py \
 -a ../../data/expression/GRCh38/expt.txt \
 -g ../../data/expression/GRCh38/Homo_sapiens.GRCh38.reference_genome_edited.fa
'''
