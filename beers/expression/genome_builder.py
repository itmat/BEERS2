from io import StringIO
import re
import os
from collections import namedtuple
from beers.expression.expression_pipeline import ExpressionPipelineException
import itertools


SingleInstanceVariant = namedtuple('SingleInstanceVariant', ['chromosome', 'position', 'description'])


class GenomeBuilderStep:

    def __init__(self, logfile, data_directory_path, parameters):
        self.data_directory_path = data_directory_path
        self.variant_line_pattern = re.compile(r'^([^|]+):(\d+) \| (.*)\tTOT')
        self.gender = parameters['gender']
        self.gender_chr_names = {name[0]: name[1] for name in zip(['X', 'Y', 'M'], parameters['gender_chr_names'])}
        self.unpaired_chr_list = self.get_unpaired_chr_list()
        self.unpaired_chr_variants = self.get_unpaired_chr_variant_data()
        self.ignore_indels = parameters['ignore_indels']
        self.ignore_snps = parameters['ignore_snps']
        self.genome_names = ['maternal', 'paternal']
        self.target_sample_name = None
        self.reference_genome = None
        self.variants_file_path = None
        self.beagle_file_path = None
        self.log_file_path = logfile

    def validate(self):
        return True

    def get_unpaired_chr_list(self):
        unpaired_chr_list = [self.gender_chr_names['M']]
        if self.gender == 'male':
            unpaired_chr_list.extend([self.gender_chr_names['X'], self.gender_chr_names['Y']])
        return unpaired_chr_list

    def get_unpaired_chr_variant_data(self):
        unpaired_chr_variants = []
        with open(self.variants_file_path) as variants_file:
            for chromosome, data in self.group_data(variants_file, lambda line: line[:line.find(':')]):
                if chromosome not in self.unpaired_chr_list:
                    continue
                for line in data:
                    match = re.match(self.variant_line_pattern, line)
                    variant_chromosome = match.group(1)
                    position = int(match.group(2)) - 1
                    variant = match.group(3).split(' | ')[0].split(":")[0]
                    unpaired_chr_variants.append(SingleInstanceVariant(variant_chromosome, position, variant))
        return unpaired_chr_variants

    def build_sequence_from_variant(self, genome, variant, reference_base):
        """
        Applies the variant provided to the custom genome provided in accordance with the variant's format (e.g.,
        D indicates delete followed by number of bases to delete, I indicates insert followed by bases to insert, and
        no D or I indicates a single base change.
        :param genome: custom genome to which the variant is applied
        :param variant: variant to apply
        :param reference_base:  base to use in place of indels when the option to ignore indels is selected.
        """

        # Indel called for but ignore indels is specified or SNP called for but ignore snps is specified,
        # then revert to reference base
        if (self.ignore_indels and ("I" in variant[0] or "D" in variant[0])) \
           or (self.ignore_snps and not ("I" in variant[0] or "D" in variant[0])):
            genome.append_segment(reference_base)
            return

        # Insert called for
        if "I" in variant[0]:
            segment_to_insert = variant[1:]
            genome.insert_segment(segment_to_insert)
            return

        # Delete called for
        if "D" in variant[0]:
            length_to_delete = int(variant[1:])
            genome.delete_segment(length_to_delete)
            return

        # SNP called for
        base_to_append = variant[0]
        genome.append_segment(base_to_append)

    def get_paired_chr_list(self):
        paired_chr_list = []
        with open(self.beagle_file_path) as beagle_file:
            for key, _ in self.group_data(beagle_file, lambda line: line.split('t')[0]):
                if key[1] != '#':
                    paired_chr_list.append(key)
        return paired_chr_list

    @staticmethod
    def group_data(lines, group_function):
        for key, values in itertools.groupby(lines, key=group_function):
            yield key, list(values)

    def locate_sample(self):
        with open(self.beagle_file_path) as beagle_file:
            for line in beagle_file:
                if line[0] == '#' and line[0:2] != "##":
                    field_headings = line.rstrip('\n').split('\t')
                    sample_names = field_headings[9:]
                    if self.target_sample_name not in sample_names:
                        raise ExpressionPipelineException(
                            f"Sample {self.target_sample_name} not found in the provided VCF file")
                    sample_index = field_headings.index(self.target_sample_name)
                    break
                elif line[0:2] == "##":
                    continue
                else:
                    raise ExpressionPipelineException(f"No sample data found.")
        return sample_index

    def execute(self, sample, reference_genome, chromosome=None):
        self.reference_genome = reference_genome
        self.chromosome_list = [chromosome] or reference_genome.keys()
        alignment_file_path = sample.input_file_path
        self.target_sample_name = sample.sample_name
        variants_filename = os.path.splitext(os.path.basename(alignment_file_path))[0] + "_variants.txt"
        self.variants_file_path = os.path.join(self.data_directory_path, variants_filename)
        variants_filename = os.path.basename(self.variants_file_path)
        self.genome_output_file_stem = \
            os.path.join(self.data_directory_path, variants_filename[:variants_filename.find('_variants')])
        self.beagle_file_path = os.path.splitext(os.path.basename(alignment_file_path))[0] + "_beagle.vcf"
        self.gender = sample.gender
        sample_index = self.locate_sample()
        paired_chr_list = self.get_paired_chr_list()
        for chromosome in self.chromosome_list:
            if chromosome in self.unpaired_chr_list:
                self.make_unpaired_chromosome(chromosome)
            elif chromosome in paired_chr_list:
                self.make_paired_chromosomes(chromosome, sample_index)
            elif not (chromosome == self.gender_chr_names['Y'] and self.gender == 'female'):
                self.make_reference_chromosome(chromosome)

    def make_reference_chromosome(self, chromosome):
        reference_sequence = self.reference_genome[chromosome]
        with open(self.log_file_path, 'a') as log_file:
            position = len(reference_sequence)
            genomes = [Genome(self.genome_names[0], chromosome, reference_sequence, position,
                              self.genome_output_file_stem),
                       Genome(self.genome_names[1], chromosome, reference_sequence, position,
                              self.genome_output_file_stem)]
            for genome in genomes:
                genome.append_segment(reference_sequence[genome.position + genome.offset:])
                log_file.write(f"Final Genome (from reference) for chromosome {chromosome}: {genome}\n")
                genome.save_to_file()

    def make_unpaired_chromosome(self, chromosome):
        reference_sequence = self.reference_genome[chromosome]
        with open(self.log_file_path, 'a') as log_file:
            genome = None
            variants = filter(lambda x: x.chromosome == chromosome, self.unpaired_chr_variants)
            for index, variant in enumerate(variants):
                if index == 0:
                    log_file.write(f"Appending"
                                   f" {len(reference_sequence[:variant.position])}"
                                   f" bases of reference sequence at reference position "
                                   f" 0 to start both genomes for chromosome {chromosome}.\n")
                    start_sequence = reference_sequence[0: variant.position]
                    if (chromosome == self.gender_chr_names['M']) or (self.gender == 'male' and chromosome ==
                                                                    self.gender_chr_names['X']):
                        genome = Genome(self.genome_names[0], chromosome, start_sequence, variant.position,
                                        self.genome_output_file_stem)
                    if chromosome == self.gender_chr_names['Y'] and self.gender == 'male':
                        genome = Genome(self.genome_names[1], chromosome, start_sequence, variant.position,
                                        self.genome_output_file_stem)
                # If the nascent genome seq position translated to reference is downstream of the variant
                # ignore the variant for this genome.
                if genome.position + genome.offset > variant.position:
                    continue
                else:
                    # If the nascent genome seq position translated to reference is upstream of the variant add
                    # the appropriate reference segment to catch up
                    if genome.position + genome.offset < variant.position:
                        log_file.write(f"Adding {variant.position - genome.position - genome.offset}"
                                       f" bases of reference sequence at reference position"
                                       f" {genome.position + genome.offset}\n")
                        genome.append_segment(reference_sequence[genome.position + genome.offset: variant.position])
                        reference_base = reference_sequence[genome.position + genome.offset: variant.position]
                        self.build_sequence_from_variant(genome, variant.description, reference_base)
            log_file.write(f"Appending"
                           f" {len(reference_sequence[genome.position + genome.offset:])}"
                           f" bases of reference sequence at reference position "
                           f" {genome.position + genome.offset} to complete the genome for"
                           f" chromosome {chromosome}.\n")
            genome.append_segment(reference_sequence[genome.position + genome.offset:])
            log_file.write(f"Final Genome for chromosome {chromosome}: {genome}\n")
            genome.save_to_file()

    def make_paired_chromosomes(self, chromosome, sample_index):
        reference_sequence = self.reference_genome[chromosome]
        with open(self.log_file_path, 'a') as log_file, open(self.beagle_file_path) as beagle_file:
            genomes = []
            for key, data in self.group_data(beagle_file, lambda line: line.split('t')):

                if key != chromosome:
                    continue

                fields = data.rstrip('\n').split('\t')
                beagle_chromosome, position, _, ref, alt, *others = fields
                sample = fields[sample_index]

                # Making position 0 based.
                position -= 1

                log_file.write(f"Appending"
                               f" {len(reference_sequence[:position])}"
                               f" bases of reference sequence at reference position "
                               f" 0 to start both genomes for chromosome {chromosome}.\n")
                start_sequence = reference_sequence[0:position]

                genomes = [Genome(self.genome_names[0], chromosome, start_sequence, position,
                                  self.genome_output_file_stem),
                           Genome(self.genome_names[1], chromosome, start_sequence, position,
                                  self.genome_output_file_stem)]

                for index, genome in enumerate(genomes):

                    # If the nascent genome seq position translated to reference is downstrem of the variant
                    # ignore the variant for this genome.
                    if genome.position + genome.offset > position:
                        continue

                    # If the nascent genome seq position translated to reference is upstream of the variant add
                    # the appropriate reference segment to catch up
                    if genome.position + genome.offset < position:
                        log_file.write(f"Adding {position - genome.position - genome.offset}"
                                       f" bases of reference sequence at reference position"
                                       f" {genome.position + genome.offset}\n")
                        genome.append_segment(reference_sequence[genome.position + genome.offset: position])

                    select_ref = sample.split("|")[index] == 0

                    if select_ref:
                        genome.append_segment(ref)
                    else:
                        if len(alt) == len(ref):
                            if self.ignore_snps:
                                genome.append_segment(ref)
                            else:
                                genome.append_segment(alt)
                        else:
                            if self.ignore_indels:
                                genome.append_segment(ref)
                            elif len(alt) > len(ref):
                                genome.insert_segment(alt)
                            else:
                                genome.delete_segment(len(ref) - len(alt))

                # Save the genome data.
                for genome in genomes:
                    log_file.write(f"Appending"
                                   f" {len(reference_sequence[genome.position + genome.offset:])}"
                                   f" bases of reference sequence at reference position "
                                   f" {genome.position + genome.offset} to complete the genome for"
                                   f" chromosome {chromosome}.\n")
                    genome.append_segment(reference_sequence[genome.position + genome.offset:])
                    log_file.write(f"Final Genome for chromosome {chromosome}: {genome}\n")
                    genome.save_to_file()


class Genome:
    """
    Holds name, chromosome, current seq, current position (0 indexed) and current offset for a nascent, custom genome.
    The current offset is such that when it is added to the current position, one arrives at the corresponding position
    (0 indexed) on the reference genome.  The object also provides methods for appending, inserting and deleting based
    upon instructions in the variants input file.
    """

    def __init__(self, name, chromosome, start_sequence, start_position, genome_output_file_stem):
        self.name = name
        self.chromosome = chromosome
        self.sequence = StringIO()
        self.sequence.write(start_sequence)
        self.position = start_position
        self.offset = 0
        self.genome_output_filename = genome_output_file_stem + '_' + self.name + ".fa"
        self.genome_indels_filename = genome_output_file_stem + '_indels_' + self.name + ".txt"
        self.indels_file = open(self.genome_indels_filename, 'a')

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
        self.indels_file.write(f"{self.chromosome}:{self.position + self.offset + 1}\tI\t{len(sequence)}\n")
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
        self.indels_file.write(f"{self.chromosome}:{self.position + self.offset + 1}\tD\t{length}\n")
        self.offset += length

    def save_to_file(self):
        """
        Saves the custom genome sequence into a single line of a fasta file.  The genome name is suffixed to the
        given output filename steam.  Since the genome sequence data is saved one chromosome at a time, the
        output file is appended to.  That means that the output file should be empty when the first chromosome
        sequence is added.  Since the sequence is memory is closed at this time, this genome can no longer be modified.
        """
        str_sequence = self.sequence.getvalue()
        self.sequence.close()
        with open(self.genome_output_filename, 'a') as genome_output_file:
            genome_output_file.write(f">{self.chromosome}\n")
            genome_output_file.write(str_sequence + "\n")

        # TODO might be a better place for this - say using a context manager?
        self.indels_file.close()

    def __str__(self):
        """
        Provide a string representation of the Genome object mainly for debugging purposes.
        :return: string representation of the Genome object
        """
        return f"name: {self.name}, chromosome: {self.chromosome}, position: {self.position}, offset: {self.offset}"
