import sys
import argparse
import os
import re
from timeit import default_timer as timer
from collections import namedtuple
from operator import attrgetter, itemgetter
import math
from io import StringIO


Variant = namedtuple('Variant', ['type', 'chromosome', 'position', 'description'])
"""
A named tuple that possesses all the attributes of a variant
type:  match (M), deletion (D), insertion (I)
chromosome: chrN
position: position on ref genome
description: description of the variant (e.g., C, IAA, D5, etc.) 
"""


class VariantsMaker:
    """
    This class creates a text file listing variants for each location in the reference genome.  The variants include
    snps and indels with the number of reads attribute to each variant.
    The text-based input file has no header and the following columns:
    1) CHROMOSOME (column 1 in a SAM file)
    2) START (column 3 in a SAM file)
    3) CIGAR  (column 4 in a SAM file)
    4) SEQ  (column 10 in a SAM file)
    The reads must have been sorted by location

    This script outputs a file that gives the full breakdown at each
    location in the genome of the number of A's, C's, G's and T's as
    well as the number of each size of insertion and deletion.
    If it's an insertion the sequence of the insertion itself is given.
    So for example a line of output like the following means
    29 reads had a C in that location, one had a T and
    also one read had an insertion TT and three reads had an insertion TTT
    chr1:10128503 | C:29 | T:1 | IT:1 ITTT:3
    """

    def __init__(self, alignment_map_filename, sort_by_entropy, depth_cutoff):
        self.entropy_sort = True if sort_by_entropy else False
        self.depth_cutoff = depth_cutoff or 10
        self.alignment_map_filename = alignment_map_filename
        self.variants_filename = os.path.splitext(alignment_map_filename)[0] + "_variants.txt"
        self.clip_at_start_pattern = re.compile("(^\d+)[SH]")
        self.clip_at_end_pattern = re.compile("\d+[SH]$")
        self.variant_pattern = re.compile("(\d+)([NMID])")
        self.indel_pattern = re.compile("\|([^|]+)")
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

    @staticmethod
    def calculate_entropy(abundances):
        """
        Use the top two abundances (if two) of the variants for the given position to compute an entropy.  If
        only one abundance is given, return 0.
        :param abundances: abundances, in terms of variant reads versus total reads of each variant at a given
        position
        :return: entropy for the given position
        """
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
        scale = 1/sum(max_abundances)
        max_abundances = [scale * max_abundance for max_abundance in max_abundances]
        return -1 * max_abundances[0] * math.log2(max_abundances[0]) - max_abundances[1] * math.log2(max_abundances[1])

    def include_annotations(self, variant_reads_per_position, total_reads_per_position, line):

        # Determine the abundances for each variant (i.e., number of reads for a variant compared with
        # the total number of reads) for this current position.
        abundances = [variant_reads_per_position[i] / total_reads_per_position
                      for i in range(0, len(variant_reads_per_position))]

        # Calculate the entropy based upon these abundances
        entropy = VariantsMaker.calculate_entropy(abundances)

        # Assemble the total reads, abundances and entropy into a string as add to the line being
        # assembled.
        annotations = f"\tTOT={total_reads_per_position}" \
                      f"\t{','.join([str(round(abundance,2)) for abundance in abundances])}" \
                      f"\tE={entropy}\n"
        line.write(annotations)
        return entropy



    def dump_to_file(self, chromosome, variant_reads, initial_write=False):
        """
        Parses the variant_reads dictionary (variant:read count) for each chromosome - position to create
        a line with the variants and their counts delimited by pipes.  Dumping each chromosome's worth of
        data at a time to avoid too sizable a dictionary.  Additionally, if the user requests a sort by entropy,
        this function will do that ordering and sent that data to stdout.
        :param chromosome: chromosome whose variants are being dumped to file
        :param variant_reads: dictionary of variants to read counts
        :param initial_write: boolean
        """
        with open(self.variants_filename, "a") as variants_file:

            # Create an in memory string to hold each line of text as it is assembled.
            line = StringIO()

            # This dictionary is only used if the user requests that the variant lines be sorted by entropy
            entropy_map = dict()

            # Record the current position being addressed.
            current_position = 0

            # Record the total number of reads found for the current position.
            total_reads_per_position = 0

            # Record the number of reads for each variant found for the current position.
            variant_reads_per_position = []

            # Iterate over the variants in the dictionary of variants to read counts sorted by the variant position
            for index, variant in enumerate(sorted(variant_reads.keys(), key=attrgetter('position'))):

                # Initial setting for the current position is the position of the first variant and we start
                # with the chromosome and position
                if index == 0:
                    current_position = variant.position
                    line.write(f'{chromosome}:{current_position}')

                # if the new position differs from the existing position, start a new line in the file and write out
                # the chromosome and new position
                if variant.position != current_position:

                    entropy = self.include_annotations(variant_reads_per_position, total_reads_per_position, line)

                    # Finally transfer the line contents to the output file
                    variants_file.write(line.getvalue())

                    # If the sort by entropy option is selected, also add to the entropy map dictionary the line
                    # entropy, keyed by the line content but only if the total number of reads exceeds the
                    # depth cutoff.
                    if self.entropy_sort and total_reads_per_position >= int(self.depth_cutoff):
                        entropy_map[line.getvalue()] = entropy

                    # Create a new in memory string to hold the next line of data
                    line = StringIO()

                    # Update the current position to the variant's position and reset the counts for total reads and
                    # reads of each variant at the current position
                    current_position = variant.position
                    total_reads_per_position = 0
                    variant_reads_per_position.clear()

                    # Write the chromosome and current position
                    line.write(f'{chromosome}:{current_position}')

                # Write the variant's description and read count and update the accumulating counters accordingly
                line.write(f' | {variant.description}:{variant_reads[variant]}')
                total_reads_per_position += variant_reads[variant]
                variant_reads_per_position.append(variant_reads[variant])

            # Handle very last position in dictionary for last chromosome
            self.include_annotations(variant_reads_per_position, total_reads_per_position, line)

            # Finally transfer the line contents to the output file
            variants_file.write(line.getvalue())

            # If the sort by entropy option is selected, also add to the entropy map dictionary the line
            # entropy, keyed by the line content but only if the total number of reads exceeds the
            # depth cutoff.
            if self.entropy_sort and total_reads_per_position >= int(self.depth_cutoff):
                entropy_map[line.getvalue()] = entropy

        # If the user selected the sort by entropy option, other the entropy_map entries in descending order
        # of entropy and print to std out.
        if self.entropy_sort:
            sorted_entropies = sorted(entropy_map.items(), key=itemgetter(1), reverse=True)
            for key, value in sorted_entropies:
                print(key, end='')

    def make_variants(self):
        """
        Iterate over the input txt file containing cigar, seq, start location, chromosome for each read and consolidate
        reads for each variant.
        """
        variant_reads = dict()
        current_chromosome = 0
        initial_write = True
        with open(self.alignment_map_filename, "r") as alignment_map_file:
            for index, line in enumerate(alignment_map_file):
                (chromosome, start, cigar, seq) = line.rstrip('\n').split('\t')
                if index == 0:
                    current_chromosome = chromosome

                # If we are starting a new chromosome, dump the data for the existing chromosome into the output
                # file and renew the variant_reads dictionary for the new chromosome.
                if current_chromosome != chromosome:
                    self.dump_to_file(current_chromosome, variant_reads, initial_write)
                    initial_write = False
                    variant_reads.clear()
                    current_chromosome = chromosome
                cigar, seq = self.remove_clips(cigar, seq)
                current_pos_in_genome = int(start)
                loc_on_read = 1

                # Iterate over the variant types and lengths in the cigar string
                for match in re.finditer(self.variant_pattern, cigar):
                    length = int(match.group(1))
                    variant_type = match.group(2)

                    # Skip over N type variants since these generally represent a read bracketed an intron
                    if variant_type == "N":
                        current_pos_in_genome += length
                        continue

                    # For a match record all the snps at the each location continuously covered by this variant type
                    if variant_type == "M":
                        stop = current_pos_in_genome + length
                        while current_pos_in_genome < stop:
                            location = current_pos_in_genome
                            variant_reads[Variant(variant_type, current_chromosome, location, seq[loc_on_read - 1])] =\
                                variant_reads.get(
                                    Variant(variant_type, current_chromosome, location, seq[loc_on_read - 1]), 0) + 1
                            loc_on_read += 1
                            current_pos_in_genome += 1
                        continue

                    # For a deletion, designate the variant named tuple description with a Dn where n is the
                    # length of the deletion starting at this position.  In this way, subsequent reads having a
                    # deletion of the same length at the same position will be added to this key.
                    if variant_type == "D":
                        location = current_pos_in_genome
                        variant_reads[Variant(variant_type, current_chromosome, location, f'D{length}')] = \
                            variant_reads.get(Variant(variant_type, current_chromosome, location, f'D{length}'), 0) + 1
                        current_pos_in_genome += length
                        continue

                    # For an insert, designate the variant name tupe description with an Ib+ where b+ are the
                    # bases to a inserted starting with this position
                    if variant_type == "I":
                        location = current_pos_in_genome
                        insertion_sequence = seq[loc_on_read - 1: loc_on_read - 1 + length]
                        variant_reads[Variant(variant_type, current_chromosome, location, f'I{insertion_sequence}')] = \
                            variant_reads.get(
                                Variant(variant_type, current_chromosome, location,  f'I{insertion_sequence}'), 0) + 1
                        loc_on_read += length

            # For each chromosome completed dump the dictionary contents to avoid a large memory footprint.
            self.dump_to_file(current_chromosome, variant_reads, initial_write)

    @staticmethod
    def main():
        """
        Entry point into the variants_maker program.  Parses the use input, created the VariantsMaker object, passing
        in the arguments and runs the process to create the variants inside a timer.
        """
        parser = argparse.ArgumentParser(description='Make Variants File')
        parser.add_argument('-a', '--alignment_map_filename',
                            help="Textfile providing chromosome, start postion, cigar,"
                                 " and sequence only for each read.")
        parser.add_argument('-s', '--sort_by_entropy', action='store_true',
                            help="Optional request to sort line in order of descreasing entropy.")
        parser.add_argument('-c', '--cutoff_depth', type=int, default=10,
                            help="Integer to indicate minimum read depth a position must have for inclusion."
                                 " If the option is not selected, a default of 10 will be applied as the minimum"
                                 " read depth.  Note that this option is used only if the sort_by_entropy option is"
                                 " invoked.")
        args = parser.parse_args()
        print(args)
        variants_maker = VariantsMaker(args.alignment_map_filename, args.sort_by_entropy, args.cutoff_depth)
        start = timer()
        variants_maker.make_variants()
        end = timer()
        sys.stderr.write(f"Variants Maker: {end - start} sec\n")


if __name__ == "__main__":
    sys.exit(VariantsMaker.main())
