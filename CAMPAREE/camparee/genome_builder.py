from io import StringIO
import re
import os
import gzip
import argparse
from collections import namedtuple
import itertools
from camparee.camparee_utils import CampareeUtils, CampareeException
from beers_utils.sample import Sample
from beers.constants import CONSTANTS


SingleInstanceVariant = namedtuple('SingleInstanceVariant', ['chromosome', 'position', 'description'])


class GenomeBuilderStep:

    def __init__(self, log_directory_path, data_directory_path, parameters=dict()):
        self.data_directory_path = data_directory_path
        self.beagle_data_file_path = os.path.join(self.data_directory_path, CONSTANTS.BEAGLE_DATA_FILE_NAME)
        self.variant_line_pattern = re.compile(r'^([^|]+):(\d+) \| (.*)\tTOT')
        self.ignore_indels = parameters.get('ignore_indels', False)
        self.ignore_snps = parameters.get('ignore_snps', False)
        self.genome_names = ['1', '2']
        self.sample_id = None
        self.reference_genome = None
        self.ploidy_data = None
        self.variants_file_path = None
        self.log_directory_path = log_directory_path

    def validate(self):
        return True

    def get_missing_chr_list(self):
        """
        Return a list of those chromosomes from chr_ploidy_data that are missing for the sample's gender.  If no
        sample gender is specified, only return a list of chromosomes from chr_ploidy_data where the chromosomes are
        missing for both genders (unlikely scenario).
        :return: list of chromosomes that are missing for this sample (likely owing to its gender)
        """
        if not self.gender:
            return [chr_ for chr_ in self.chr_ploidy_data.keys()
                    if self.chr_ploidy_data[chr_][CONSTANTS.MALE_GENDER] !=
                    self.chr_ploidy_data[chr_][CONSTANTS.FEMALE_GENDER]]
        return [chr_ for chr_ in self.chr_ploidy_data.keys() if self.chr_ploidy_data[chr_][self.gender] == 0]

    def get_unpaired_chr_list(self):
        """
        Return a list of those chromosomes from chr_ploidy_data that are unpaired for the sample's gender.  If no
        sample gender is specified, only return a list of chromosomes from chr_ploidy_data where the chromosomes are
        unpaired for both genders.
        :return: list of chromosomes that are unpaired for this sample (likely owing to its gender)
        """
        if not self.gender:
            return [chr_ for chr_ in self.chr_ploidy_data.keys()
                    if self.chr_ploidy_data[chr_][CONSTANTS.MALE_GENDER] ==
                    self.chr_ploidy_data[chr_][CONSTANTS.FEMALE_GENDER] == 1]
        return [chr_ for chr_ in self.chr_ploidy_data.keys() if self.chr_ploidy_data[chr_][self.gender] == 1]

    def get_paired_chr_list(self):
        """
        Return a list of those chromosomes from chr_ploidy_data that are paired for the sample's gender.  If no
        sample gender is specified, only return a list of chromosomes from chr_ploidy_data where the chromosomes are
        paired for both genders.
        :return: list of chromosomes that are paired for this sample (likely owing to its gender)
        """
        if not self.gender:
            return [chr_ for chr_ in self.chr_ploidy_data.keys()
                    if self.chr_ploidy_data[chr_][CONSTANTS.MALE_GENDER] ==
                    self.chr_ploidy_data[chr_][CONSTANTS.FEMALE_GENDER] == 2]
        return [chr_ for chr_ in self.chr_ploidy_data.keys() if self.chr_ploidy_data[chr_][self.gender] == 2]

    def get_unpaired_chr_variant_data(self):
        """
        There should be at most, one variant for any given position in an unpaired chromosome.  This method groups
        the variant records by chromosome for those chromosomes found in the unpaired chr list and adds a single
        instance variant to the an unpaired_chr_variants list for every such variant found and returns the list.
        :return: A list of all unpaired chromosome variants
        """
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

    @staticmethod
    def group_data(lines, group_function):
        """
        Returns data grouped by the provided function
        :param lines:  the lines of data to be grouped
        :param group_function: The function to apply to determine the
        groupping.
        :return: a generator providing the next key (the groupping
        parameter) and the groupped data as a list.
        """
        for key, values in itertools.groupby(lines, key=group_function):
            yield key, list(values)

    def locate_sample(self):
        """
        Find the position of the sample in the beagle data
        :return: The position of the sample in a line of beagle data
        """
        with gzip.open(self.beagle_data_file_path) as beagle_file:
            for line in beagle_file:
                line = line.decode('ascii')
                if line[0] == '#' and line[0:2] != "##":
                    field_headings = line.rstrip('\n').split('\t')
                    sample_ids = field_headings[9:]
                    if self.sample_id not in sample_ids:
                        raise CampareeException(f"Sample {self.sample_id} not found in the provided VCF file")
                    sample_index = field_headings.index(self.sample_id)
                    break
                elif line[0:2] == "##":
                    continue
                else:
                    raise CampareeException(f"No sample data found.")
        return sample_index

    def execute(self, sample, chr_ploidy_data, reference_genome, chromosome_list=[]):
        """
        Entry point for genome builder.  Uses the chr_ploidy_data and reference_genome resources along with beagle and
        variant finder output to build two custom genomes.
        :param sample: sample for which the genome is being built
        :param chr_ploidy_data: dictionary indicating chromosomes to be processed and their ploidy based on sample
        gender.
        :param reference_genome: dictionary relating chr to its reference sequence
        :param chromosome_list: a debug feature that overrides the chr_ploidy_data chr list.  Useful for testing a
        specific chromosome only or a small subset of chromosomes.
        """
        self.chr_ploidy_data = chr_ploidy_data
        self.reference_genome = reference_genome
        # The chromosome list derived from the chr_ploidy_data is the gold standard.  Only those chromosomes/contigs
        # are processed.
        self.chromosome_list = chromosome_list or list(chr_ploidy_data.keys())
        self.sample_id = f'sample{sample.sample_id}'
        self.gender = sample.gender
        self.variants_file_path = os.path.join(self.data_directory_path, self.sample_id, CONSTANTS.VARIANTS_FILE_NAME)
        sample_index = self.locate_sample()
        self.genome_output_file_stem = os.path.join(self.data_directory_path, self.sample_id, 'custom_genome')
        self.log_file_path = os.path.join(self.log_directory_path, self.sample_id, __class__.__name__ + ".log")
        self.unpaired_chr_list = self.get_unpaired_chr_list()
        self.unpaired_chr_variants = self.get_unpaired_chr_variant_data()
        self.paired_chr_list = self.get_paired_chr_list()
        # How the chromosome should be processed depends upon whether it is paired or not.
        for chromosome in self.chromosome_list:
            if chromosome in self.unpaired_chr_list:
                self.make_unpaired_chromosome(chromosome)
            elif chromosome in self.paired_chr_list:
                self.make_paired_chromosome(chromosome, sample_index)
            elif chromosome not in self.get_missing_chr_list() and chromosome in self.chr_ploidy_data.keys():
                # We should never get here since the conditions above should be mutually exclusive.
                self.make_reference_chromosome(chromosome)

    def make_reference_chromosome(self, chromosome):
        """
        Here, the reference sequence for the given chromosome is copied as it, into the custom
        genomes.
        :param chromosome: The chromosome for which the reference sequence is used.
        """
        print(f'Processing chromosome {chromosome} using reference genome chromosome sequence.')
        reference_sequence = self.reference_genome[chromosome]
        with open(self.log_file_path, 'a') as log_file:
            position = len(reference_sequence)
            genomes = [Genome(self.genome_names[0], chromosome, reference_sequence, position,
                              self.genome_output_file_stem)]

            # Add reference sequence to both genomes only if the chromosome is designated as paired.
            if chromosome in self.paired_chr_list:
                genomes.append(Genome(self.genome_names[1], chromosome, reference_sequence, position,
                                      self.genome_output_file_stem))
            for genome in genomes:
                log_file.write(f"Final Genome (from reference) for chromosome {chromosome}: {genome}\n")
                genome.save_to_file()

    def make_unpaired_chromosome(self, chromosome):
        """
        Here, the samples variants data is threaded together with the reference sequence to create a custom
        sequence for the given chromosome.
        :param chromosome: The chromosome for which the reference sequence is altered by variant data.
        """
        reference_sequence = self.reference_genome[chromosome]
        genome = None
        variants = list(filter(lambda x: x.chromosome == chromosome, self.unpaired_chr_variants))
        # If no variants are found for this chromosome, copy over the reference sequence instead.
        if not variants:
            self.make_reference_chromosome(chromosome)
            return
        print(f'Processing chromosome {chromosome} from unpaired chromosome list.')
        with open(self.log_file_path, 'a') as log_file:
            for index, variant in enumerate(variants):
                if index == 0:
                    log_file.write(f"Appending"
                                   f" {len(reference_sequence[:variant.position])}"
                                   f" bases of reference sequence at reference position "
                                   f" 0 to start both genomes for chromosome {chromosome}.\n")
                    start_sequence = reference_sequence[0: variant.position]

                    # We don't care where an unpaired chromosome ends up.  So we put it into the first of the two
                    # genomes.  Whether the chromosome was contributed by the mother or the father is of no
                    # importance.
                    genome = Genome(self.genome_names[0], chromosome, start_sequence, variant.position,
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
                                       f" {genome.position + genome.offset} for chromosome {chromosome}\n")
                        genome.append_segment(reference_sequence[genome.position + genome.offset: variant.position])
                        reference_base = reference_sequence[genome.position + genome.offset: variant.position]
                        self.build_sequence_from_variant(genome, variant.description, reference_base)
                        log_file.write(f"Currently: {genome}\n")
            log_file.write(f"Appending"
                           f" {len(reference_sequence[genome.position + genome.offset:])}"
                           f" bases of reference sequence at reference position "
                           f" {genome.position + genome.offset} to complete the genome for"
                           f" chromosome {chromosome}.\n")
            genome.append_segment(reference_sequence[genome.position + genome.offset:])
            log_file.write(f"Final Genome for chromosome {chromosome}: {genome}\n")
            genome.save_to_file()

    def make_paired_chromosome(self, chromosome, sample_index):
        """
        Here, the beagle data for the given sample is threaded together with the reference sequence to create a
        custom sequence for the given chromosome
        :param chromosome: The chromosome for which the reference sequence is altered by beagle data.
        :param sample_index: identifies the position of the subject sample in the beagle data.
        """
        reference_sequence = self.reference_genome[chromosome]
        with gzip.open(self.beagle_data_file_path) as beagle_file:
            for key, data in self.group_data(beagle_file, lambda line: line.decode('ascii').split('\t')[0]):
                if key == chromosome:
                    break
            else:
                # If the paired chromosome is not found in the beagle file, use the reference chromosome.
                self.make_reference_chromosome(chromosome)
                return

            print(f'Processing chromosome {chromosome} from paired chromosome list.')
            with open(self.log_file_path, 'a') as log_file:

                genomes = [Genome(self.genome_names[0], chromosome, '', 0, self.genome_output_file_stem),
                           Genome(self.genome_names[1], chromosome, '', 0, self.genome_output_file_stem)]

                for datum in data:

                    fields = datum.decode('ascii').rstrip('\n').split('\t')
                    log_file.write(f"Beagle Field: {fields}\n")
                    beagle_chromosome, position, _, ref, alts, *others = fields
                    position = int(position)
                    alts = alts.split(',')
                    sample = fields[sample_index].strip().split('|')
                    alt_indexes = [int(s) for s in sample]

                    # Making position 0 based.
                    position -= 1

                    for genome_index, genome in enumerate(genomes):

                        # If the nascent genome seq position translated to reference is downstream of the variant
                        # ignore the variant for this genome.
                        if genome.position + genome.offset > position:
                            continue

                        # If the nascent genome seq position translated to reference is upstream of the variant add
                        # the appropriate reference segment to catch up
                        if genome.position + genome.offset < position:
                            log_file.write(f"Adding {position - genome.position - genome.offset}"
                                           f" bases of reference sequence at reference position"
                                           f" {genome.position + genome.offset} to genome_{genome.name}\n")
                            genome.append_segment(reference_sequence[genome.position + genome.offset: position])

                        alt_index = alt_indexes[genome_index]
                        if alt_index == 0:
                            genome.append_segment(ref)
                        else:
                            alt = alts[alt_index-1]

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
                                    # TODO: does this always correctly "delete" the right segment?
                                    #  Eg: if ref is ACGT and alt is CG, then there'd actually be two deletions
                                    #  Unclear if we actually get these, but also whether deletions are at the start or end
                                    genome.delete_segment(len(ref) - len(alt))
                            log_file.write(f"Currently: {genome}\n")

                # Save the genome data.
                for genome in genomes:
                    log_file.write(f"Appending"
                                   f" {len(reference_sequence[genome.position + genome.offset:])}"
                                   f" bases of reference sequence at reference position "
                                   f" {genome.position + genome.offset} to complete genome_{genome.name} for"
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


if __name__ == "__main__":
    import pandas as pd
    parser = argparse.ArgumentParser(description='Make Genome Files')
    parser.add_argument('-r', '--reference_genome_file_path',
                        help="Fasta file containing the reference genome")
    parser.add_argument('-l', '--log_directory_path', help="Path to log directory.")
    parser.add_argument('-x', '--gender', action='store', choices=['male', 'female'],
                        help="Gender of input sample (male or female).")
    parser.add_argument('-c', '--chromosomes', type=lambda chrs: [chr_ for chr_ in chrs.split(',')],
                        help="optional chromosome list")
    parser.add_argument('-i', '--ignore_indels', action='store_true',
                        help="Use the reference genome base in place of an indel.  Defaults to false.")
    parser.add_argument('-j', '--ignore_snps', action='store_true',
                        help="Use the reference genome base in place of a snp.  Defaults to False.")
    parser.add_argument('-d', '--data_directory_path', help='Path to data directory')
    parser.add_argument('-s', '--sample_id', type=int, help='sample name in vcf when prepended with sample')
    parser.add_argument('-p', '--chr_ploidy_file_path', type=str, help="Path to chromosome ploidy file")
    args = parser.parse_args()
    print(args)

    test_parameters = {"ignore_snps": args.ignore_snps,
                       "ignore_indels": args.ignore_indels}
    test_sample = Sample(args.sample_id, "debug sample", None, None, args.gender)
    reference_genome_ = CampareeUtils.create_genome(args.reference_genome_file_path)
    chr_ploidy_data_ = CampareeUtils.create_chr_ploidy_data(args.chr_ploidy_file_path)

    # Remove old genome data and log files if present
    sample_data_folder = os.path.join(args.data_directory_path, f'sample{test_sample.sample_id}')
    for item in os.listdir(sample_data_folder):
        if item.startswith("custom_genome"):
            os.remove(os.path.join(sample_data_folder, item))
    sample_log_folder = os.path.join(args.log_directory_path, f'sample{test_sample.sample_id}')
    for item in os.listdir(sample_log_folder):
        if item.startswith("GenomeBuilderStep"):
            os.remove(os.path.join(sample_log_folder, item))

    genome_builder = GenomeBuilderStep(args.log_directory_path, args.data_directory_path, test_parameters)
    genome_builder.execute(test_sample, chr_ploidy_data_, reference_genome_, args.chromosomes)

    results = CampareeUtils.compare_genome_sequence_lengths(args.reference_genome_file_path,
                                                            os.path.join(sample_data_folder, 'custom_genome_1.fa'),
                                                            os.path.join(sample_data_folder, 'custom_genome_2.fa'),
                                                            args.chromosomes or chr_ploidy_data_.keys())
    df = pd.DataFrame.from_dict(results, orient='index', columns=['Reference Genome', 'Genome 1', 'Genome 2'])
    print(df)

'''
Example
python genome_builder.py \
-c "MT" \
-s 1 \
-d ../../data/pipeline_results_run99/expression_pipeline/data \
-r ../../resources/index_files/GRCh38/Homo_sapiens.GRCh38.reference_genome.fa.gz \
-p ../../resources/index_files/GRCh38/Homo_sapiens.GRCh38.chr_ploidy.txt \
-l ../../data/pipeline_results_run99/expression_pipeline/logs \
-x female

python genome_builder.py \
-s 1 \
-d ../../data/pipeline_results_run89/expression_pipeline/data \
-r ../../resources/index_files/GRCh38/Homo_sapiens.GRCh38.reference_genome.fa.gz \
-p ../../resources/index_files/GRCh38/Homo_sapiens.GRCh38.chr_ploidy.txt \
-l ../../data/pipeline_results_run89/expression_pipeline/logs \
-x female
'''
