import sys
import argparse
import os
import re
from timeit import default_timer as timer
from collections import namedtuple
from operator import attrgetter


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

    def __init__(self, alignment_map_filename):
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
        clip_at_start = re.search(self.clip_at_start_pattern, cigar)
        if clip_at_start:
            cigar = re.sub(self.clip_at_start_pattern, "", cigar)
            seq = seq[int(clip_at_start.group(1)):]
        cigar = re.sub(self.clip_at_end_pattern, "", cigar)
        return cigar, seq

    def dump_to_file(self, chromosome, output):
        with open(self.variants_filename, "a") as variants_file:
            current_position = 0
            for index, variant in enumerate(sorted(output.keys(), key=attrgetter('position'))):
                if index == 0:
                    current_position = variant.position
                if variant.position != current_position:
                    variants_file.write("\n")
                    current_position = variant.position
                    variants_file.write(f'{chromosome}:{current_position}')
                variants_file.write(f' | {variant.description}:{output[variant]}')

    def make_variants(self):
        output = dict()
        current_chromosome = 0
        with open(self.alignment_map_filename, "r") as alignment_map_file:
            for index, line in enumerate(alignment_map_file):
                (chromosome, start, cigar, seq) = line.rstrip('\n').split('\t')
                if index == 0:
                    current_chromosome = chromosome
                if current_chromosome != chromosome:
                    self.dump_to_file(chromosome, output)
                    output.clear()
                    current_chromosome = chromosome
                cigar, seq = self.remove_clips(cigar, seq)
                current_pos_in_genome = int(start)
                loc_on_read = 1
                for match in re.finditer(self.variant_pattern, cigar):
                    length = int(match.group(1))
                    variant_type = match.group(2)
                    if variant_type == "N":
                        current_pos_in_genome += length
                        continue
                    if variant_type == "M":
                        stop = current_pos_in_genome + length
                        while current_pos_in_genome < stop:
                            location = current_pos_in_genome
                            output[Variant(variant_type, current_chromosome, location, seq[loc_on_read - 1])] =\
                                output.get(
                                    Variant(variant_type, current_chromosome, location, seq[loc_on_read - 1]), 0) + 1
                            loc_on_read += 1
                            current_pos_in_genome += 1
                        continue
                    if variant_type == "D":
                        location = current_pos_in_genome
                        output[Variant(variant_type, current_chromosome, location, f'D{length}')] = \
                            output.get(Variant(variant_type, current_chromosome, location, f'D{length}'), 0) + 1
                        current_pos_in_genome += length
                        continue
                    if variant_type == "I":
                        location = current_pos_in_genome
                        insertion_sequence = seq[loc_on_read - 1: loc_on_read - 1 + length]
                        output[Variant(variant_type, current_chromosome, location, f'I{insertion_sequence}')] = \
                            output.get(
                                Variant(variant_type, current_chromosome, location,  f'I{insertion_sequence}'), 0) + 1
                        loc_on_read += length
            self.dump_to_file(current_chromosome, output)

    @staticmethod
    def main():
        parser = argparse.ArgumentParser(description='Make Variants File')
        parser.add_argument('-a', '--alignment_map_filename')
        args = parser.parse_args()
        variants_maker = VariantsMaker(args.alignment_map_filename)
        start = timer()
        variants_maker.make_variants()
        end = timer()
        print(f"Variants Maker: {end - start}")


if __name__ == "__main__":
    sys.exit(VariantsMaker.main())
