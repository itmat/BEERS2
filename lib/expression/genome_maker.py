from io import StringIO
import numpy as np
import sys
import re
import os
import argparse
from timeit import default_timer as timer
from operator import itemgetter


class GenomeMaker:

    def __init__(self, variants_filename, genome_ref_filename, genome_output_file_stem, seed, threshold, log_filename):
        self.variants_filename = variants_filename
        self.genome_output_file_stem = genome_output_file_stem
        self.genome_ref_filename = genome_ref_filename
        self.log_filename = log_filename
        np.random.seed(seed)
        self.abundance_threshold = threshold
        self.variant_line_pattern = re.compile('^([^|]+):(\d+) \| (.*)\tTOT')
        self.genome_names = ['seq1', 'seq2']
        for genome_name in self.genome_names:
            try:
                os.remove(self.genome_output_file_stem + "_" + genome_name + ".fa")
            except OSError:
                pass

    def get_most_abundant_variants(self, variants):
        max_variants = []
        variant_reads = {variant.split(':')[0]: int(variant.split(':')[1]) for variant in variants}
        if len(variant_reads) == 1:
            return [next(iter(variant_reads.keys()))]
        for _ in range(2):
            max_variant = max(variant_reads.items(), key=itemgetter(1))[0]
            max_variants.append((max_variant, variant_reads[max_variant]))
            del variant_reads[max_variant]
        total_reads = sum(read for _, read in max_variants)
        if max_variants[0][1]/total_reads < self.abundance_threshold:
            return [max_variants[0][0]]
        return [variant for variant, read in max_variants]

    @staticmethod
    def build_sequence_from_variant(genome, variant):

        # Insert called for
        if "I" in variant[0]:
            segment_to_insert = variant[1:]
            genome.insert_segment(segment_to_insert)

        # Delete called for
        elif "D" in variant[0]:
            length_to_delete = int(variant[1:])
            genome.delete_segment(length_to_delete)

        # SNP called for
        else:
            base_to_append = variant[0]
            genome.append_segment(base_to_append)

    def make_genome(self):

        with open(self.log_filename, 'w') as log_file:
            with open(self.variants_filename) as variants_file, open(self.genome_ref_filename) as genome_ref_file:
                reference_chromosome = None
                reference_sequence = None
                building_chromosome = False
                genomes = list()
                for variant_line in variants_file:
                    match = re.match(self.variant_line_pattern, variant_line)
                    chromosome = match.group(1)

                    # Substraction to account for fact that variant text file indexed at 1
                    variant_position = int(match.group(2)) - 1
                    variants = match.group(3).split(' | ')

                    # Starting new chromosome (possibly the first one)
                    while not reference_chromosome or reference_chromosome < chromosome:

                        # At least one chromosome already completed - finish and save
                        if building_chromosome:
                            for genome in genomes:
                                log_file.write(f"Appending"
                                               f" {len(reference_sequence[genome.position + genome.offset:])}\n")
                                genome.append_segment(reference_sequence[genome.position + genome.offset:])
                                log_file.write(f"{genome}\n")
                                genome.save_to_file(self.genome_output_file_stem)

                        # Get a reference sequence for the next chromosome in the reference genome file
                        line = genome_ref_file.readline()
                        if line.startswith(">"):
                            reference_chromosome = line[1:].rstrip()

                            # Only set up genomes if the reference chromosome is the same at the one currently
                            # being read from the variants file.  Otherwise, skip over and get the next chromosome from
                            # the reference genome file.
                            if reference_chromosome == chromosome:
                                building_chromosome = True
                                reference_sequence = genome_ref_file.readline().rstrip()
                                genomes = list()
                                start_sequence = reference_sequence[0: variant_position]
                                genomes.append(
                                    Genome(self.genome_names[0], chromosome, start_sequence, variant_position))
                                genomes.append(
                                    Genome(self.genome_names[1], chromosome, start_sequence, variant_position))
                            else:
                                genome_ref_file.readline()
                                building_chromosome = False

                    # Return the top two (based on number of reads) of the variants on this line.
                    max_variants = self.get_most_abundant_variants(variants)

                    if max_variants[0] != reference_sequence[variant_position]:
                        log_file.write(f"Reference at {chromosome}:{variant_position} is "
                                       f"{reference_sequence[variant_position]} whereas max variant "
                                       f"is {max_variants[0]}\n")

                    for genome in genomes:

                        # If the nascent genome seq position translated to reference is downstrem of the variant
                        # ignore the variant for this genome.
                        if genome.position + genome.offset > variant_position:
                            continue
                        else:
                            # If the nascent genome seq position translated to reference is upstream of the variant add
                            # the appropriate reference segment to catch up
                            if genome.position + genome.offset < variant_position:
                                log_file.write(f"Adding {variant_position - genome.position - genome.offset}"
                                               f" bases of reference sequence at reference position"
                                               f" {genome.position + genome.offset}\n")
                                genome.append_segment(
                                    reference_sequence[genome.position + genome.offset: variant_position])

                            # If only one variant, apply directly to genome
                            if len(max_variants) == 1:
                                self.build_sequence_from_variant(genome, max_variants[0])

                            # If two variants exist, toss a coin for one
                            else:
                                max_variant = np.random.choice([max_variants[0], max_variants[1]], p=[0.5, 0.5])
                                self.build_sequence_from_variant(genome, max_variant)
                                max_variants.remove(max_variant)

                # Variants file exhausted, must be done with last chromosome.
                for genome in genomes:
                    genome.append_segment(reference_sequence[genome.position + genome.offset:])
                    log_file.write(f"{genome}\n")
                    genome.save_to_file(self.genome_output_file_stem)

    @staticmethod
    def main():
        parser = argparse.ArgumentParser(description='Make Genome Files')
        parser.add_argument('-v', '--variants_filename',
                            help="Textfile providing variants as provided by the ouput"
                                 " of the variants_maker.py script.")
        parser.add_argument('-r', '--reference_genome_filename',
                            help="Fasta file containing the reference genome")
        parser.add_argument('-g', '--genome_output_file_stem',
                            help="Path of output genome file.  Will be suffixed _seq1.fa or _seq2.fa")
        parser.add_argument('-s', '--seed', type=int, default=100,
                            help="Integer to be used as a seed for the random number generator."
                                 "  Value defaults to 100.")
        parser.add_argument('-t', '--threshold', type=float, default=0.03,
                            help="Abundance threshold of alt allele.  Defaults to 0.03")
        parser.add_argument('-l', '--log_filename', help="Log file.")
        args = parser.parse_args()
        print(args)
        genome_maker = GenomeMaker(args.variants_filename,
                                   args.reference_genome_filename,
                                   args.genome_output_file_stem,
                                   args.seed,
                                   args.threshold,
                                   args.log_filename)
        start = timer()
        genome_maker.make_genome()
        end = timer()
        sys.stderr.write(f"Genome Maker: {end - start} sec\n")


class Genome:
    """
    Holds name, chromosome, current sequence, current position (0 indexed) and current offset for a nascent, custom genome.
    The current offset is such that when it is added to the current position, one arrives at the corresponding position
    (0 indexed) on the reference genome.  The object also provides methods for appending, inserting and deleting based
    upon instructions in the variants input file.
    """

    def __init__(self, name, chromosome, start_sequence, start_position):
        self.name = name
        self.chromosome = chromosome
        self.sequence = StringIO()
        self.sequence.write(start_sequence)
        self.position = start_position
        self.offset = 0

    def append_segment(self, sequence):
        """
        Append the given sequence segment to the custom genome.  Since the sequence segment either has a one to one
        correspondence with that of reference genome or is a sequence segment drawn from the reference genome;
        execution of this method does not alter the current position of the custom genome relative to the current
        position of the reference genome/variant.  So position advances by the sequence segment length but offset
        remains unchangeed.
        :param sequence: sequence segment to append
        """
        self.sequence.write(sequence)
        self.position += len(sequence)

    def insert_segment(self, sequence):
        """
        Insert the given sequence segment into the custom genome.  Since the given sequence segment does not correspond
        to anything in the reference genome; the current position of the custom genome relative to the current position
        of the reference genome/variant does change by the length of the sequence segment.  Since the custom genome
        sequence is advancing while the reference sequence is not, the sequence segment length is subtracted from the
        offset while the genome current position is advanced by the length of the sequence segment.
        :param sequence: sequence segment to insert
        """
        self.sequence.write(sequence)
        self.position += len(sequence)
        self.offset += -1 * len(sequence)

    def delete_segment(self, length):
        """
        Skip over (delete) a length of the reference sequence.  Since the reference sequence is advancing while the
        custom sequence is not, the relative current position of the genome again changes relative to the current
        position of the reference sequence.  As such, the current genome position does not advance but the offset
        increases by the length provided.
        :param length: number of bases in the reference sequence to skip over.
        """
        self.offset += length

    def save_to_file(self, genome_output_file_stem):
        str_sequence = self.sequence.getvalue()
        self.sequence.close()
        with open(genome_output_file_stem + '_' + self.name + ".fa", 'a') as genome_file:
            genome_file.write(f">{self.chromosome}\n")
            genome_file.write(str_sequence + "\n")

    def __str__(self):
        return f"name: {self.name}, chromosome: {self.chromosome}, position: {self.position}, offset: {self.offset}"


if __name__ == "__main__":
    sys.exit(GenomeMaker.main())


'''
Example:

python genome_maker.py -v  ../../data/preBEERS/ETAM080_grp1.gene.norm.chr21_22.all_unique_mappers.fw_only_variants.txt \
-g ../../data/preBEERS/human_21_22_genome -r ../../data/preBEERS/hg19_chr21_22_ref_edited.fa \
-l ../../data/preBEERS/genome_maker.log
'''
