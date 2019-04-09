from io import StringIO
import re
import os
import sys
import gzip
import argparse
from collections import namedtuple
import itertools
from camparee.camparee_utils import CampareeUtils, CampareeException
from camparee.abstract_camparee_step import AbstractCampareeStep
from beers_utils.sample import Sample
from beers_utils.constants import CONSTANTS


SingleInstanceVariant = namedtuple('SingleInstanceVariant', ['chromosome', 'position', 'description'])


class GenomeBuilderStep(AbstractCampareeStep):

    #Prefix for names of all (non-log) files output by this script.
    GENOME_BUILDER_OUTPUT_PREFIX = "custom_genome"
    #Name of file where script logging is stored
    GENOME_BUILDER_LOG_FILENAME = "GenomeBuilderStep.log"

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

    def execute(self, sample, chr_ploidy_data, reference_genome, chromosome_list=None):
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
        self.chromosome_list = chromosome_list if chromosome_list else list(chr_ploidy_data.keys())
        self.sample_id = f'sample{sample.sample_id}'
        self.gender = sample.gender
        self.variants_file_path = os.path.join(self.data_directory_path, self.sample_id,
                                               CONSTANTS.VARIANTS_FILE_NAME)
        sample_index = self.locate_sample()
        self.genome_output_file_stem = os.path.join(self.data_directory_path,
                                                    self.sample_id,
                                                    GenomeBuilderStep.GENOME_BUILDER_OUTPUT_PREFIX)
        self.log_file_path = os.path.join(self.log_directory_path, self.sample_id,
                                          GenomeBuilderStep.GENOME_BUILDER_LOG_FILENAME)
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
        #Final entry in log file to indicate execution() method finished.
        with open(self.log_file_path, 'a') as log_file:
            log_file.write('\nALL DONE!\n')

    def get_commandline_call(self, sample, chr_ploidy_file_path, reference_genome_file_path, chromosome_list=None):
        """
        Prepare command to execute the GenomeBuilderStep from the command line,
        given all of the arugments used to run the execute() function.

        Parameters
        ----------
        sample : Sample
            Sample for which to construct parental genomes
        chr_ploidy_file_path : string
            File that maps chromosome names to their male/female ploidy.
        reference_genome_file_path : string
            File that maps chromosome names in reference to nucleotide sequence.
        chromosome_list : list
            A debug feature that overrides the chr_ploidy_data chr list. Useful
            for testing a specific chromosome only or a small subset of
            chromosomes.

        Returns
        -------
        string
            Command to execute on the command line. It will perform the same
            operations as a call to execute() with the same parameters.

        """
        #Retrieve path to the variants_finder.py script.
        genome_builder_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        genome_builder_path = genome_builder_path.rstrip('c')

        #Note, the construction of this command line call deviates a bit from the
        #other steps (various options are specified explicitly, rather than passing
        #a json dump of the config parameters). This command was already in place
        #before we added command line functionality to the majority of the other
        #steps. I'm leaving this as-is to maintain functionalty with anyone else's
        #previous code or tests.

        command = (f" python {genome_builder_path}"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --ignore_indels {self.ignore_indels}"
                   f" --ignore_snps {self.ignore_snps}"
                   f" --sample '{repr(sample)}'"
                   f" --chr_ploidy_file_path {chr_ploidy_file_path}"
                   f" --reference_genome_file_path {reference_genome_file_path}")

        if chromosome_list:
            #Create comma-searated list of chromosomes
            command += f" --chromosomes {','.join(str(chr) for chr in chromosome_list)}"

        return command

    def get_validation_attributes(self, sample):
        """
        Prepare attributes required by is_output_valid() function to validate
        output generated the GenomeBuilderStep job corresponding to the given
        sample.

        Parameters
        ----------
        sample : Sample
            Sample for which custom parental genomes will be generated.

        Returns
        -------
        dict
            A GenomeBuilderStep job's data_directory, log_directory, sample_id,
            and a list of the genome names used by the GenomeBuilderStep to refer
            to each of the parental genomes (i.e. 1 and 2 for male and female
            parent, respectively).
        """
        validation_attributes = {}
        validation_attributes['data_directory'] = self.data_directory_path
        validation_attributes['log_directory'] = self.log_directory_path
        validation_attributes['sample_id'] = sample.sample_id
        validation_attributes['genome_names'] = self.genome_names
        return validation_attributes

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
                    # At this point, the nascent genome seq position is at that of the variant.  So change introduced
                    # by the variant can be incorporated.
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
        custom sequence for the given chromosome.  Below is a snippet of a beagle vcf file for 6 samples along with
        a header.

        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1 sample2 sample3 sample4 sample5 sample6
        ...
        chr1    257558  .       A       G       .       PASS    .       GT      0|1     0|0     0|0     0|0     0|0     0|0
        chr1    257559  .       G       C,GAG   .       PASS    .       GT      0|1     0|2     0|0     0|0     0|0     0|0
        chr1    257560  .       C       A       .       PASS    .       GT      1|0     1|0     0|0     1|0     0|0     0|0
        chr1    257570  .       C       CAA,CA  .       PASS    .       GT      1|0     0|0     0|0     0|0     0|1     2|0

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

                    # Collect alt, ref selections for this sample where 0 = ref and > 0 = alt (i.e., 0/1 = ref for 1st
                    # parent chr and alt for 2nd while 1/2 = 1st alt for 1st parent chr and 2nd alt for 2nd, etc.)
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

                        # Identify alt or ref for current parent
                        alt_index = alt_indexes[genome_index]

                        # Apply the ref or alt as appropriate
                        if alt_index == 0:
                            genome.append_segment(ref)
                        else:
                            # Get the appropriate alt as given by the index - 1 (since the 1st alt is given by 1)
                            alt = alts[alt_index-1]

                            # Apply a snp alt only if the ignore_snps parameter is not set
                            if len(alt) == len(ref):
                                if self.ignore_snps:
                                    genome.append_segment(ref)
                                else:
                                    genome.append_segment(alt)
                            # Otherwise, apply only if the ignore_indels parameter is not set
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

    @staticmethod
    def is_output_valid(validation_attributes):
        """
        Check if output of GenomeBuilderStep for a specific job/execution is
        correctly formed and valid, given a job's data directory, log directory,
        and sample id. Prepare these attributes for a given sample's jobs using
        the get_validation_attributes() method.

        Parameters
        ----------
        validation_attributes : dict
            A job's data_directory, log_directory, sample_id, and the list of
            genome names used by the GenomeBuilderStep to refer to each of the
            parental genomes (i.e. 1 and 2 for male and female parent, respectively).

        Returns
        -------
        boolean
            True  - GenomeBuilderStep output files were created and are well formed.
            False - GenomeBuilderStep output files do not exist or are missing data.

        """
        data_directory = validation_attributes['data_directory']
        log_directory = validation_attributes['log_directory']
        sample_id = validation_attributes['sample_id']
        genome_names = validation_attributes['genome_names']

        valid_output = False

        genome_builder_outfile_path = os.path.join(data_directory, f"sample{sample_id}",
                                                   GenomeBuilderStep.GENOME_BUILDER_OUTPUT_PREFIX)
        genome_builder_logfile_path = os.path.join(log_directory, f"sample{sample_id}",
                                                   GenomeBuilderStep.GENOME_BUILDER_LOG_FILENAME)

        if os.path.isfile(genome_builder_logfile_path):
            #Read last line in genome_builder log file
            line = ""
            with open(genome_builder_logfile_path, "r") as genome_builder_log_file:
                for line in genome_builder_log_file:
                    line = line.rstrip()
            if line == "ALL DONE!":
                #Now check all output files were created for each of the custom
                #parental genomes. There are multiple files per genome. The naming
                #schemes for these output files are defined in the Genome class
                #below.
                custom_genome_output_list = []

                for name in genome_names:
                    genome_output_filename = \
                        Genome.GENOME_OUTPUT_FILENAME_PATTERN.format(
                            genome_output_file_stem=genome_builder_outfile_path,
                            genome_name=name
                        )
                    indel_output_filename = \
                        Genome.INDEL_OUTPUT_FILENAME_PATTERN.format(
                            genome_output_file_stem=genome_builder_outfile_path,
                            genome_name=name
                        )
                    genome_mapper_filename = \
                        Genome.GENOME_MAPPER_FILENAME_PATTERN.format(
                            genome_output_file_stem=genome_builder_outfile_path,
                            genome_name=name
                        )
                    custom_genome_output_list.append(os.path.isfile(genome_output_filename))
                    custom_genome_output_list.append(os.path.isfile(indel_output_filename))
                    custom_genome_output_list.append(os.path.isfile(genome_mapper_filename))

                if all(custom_genome_output_list):
                    valid_output = True

        return valid_output

    @staticmethod
    def main():
        """
        Entry point into script. Allows script to be executed/submitted via the
        command line.
        """

        parser = argparse.ArgumentParser(description='Make Genome Files')

        parser.add_argument('-l', '--log_directory_path', required=True,
                            help="Path to log directory.")
        parser.add_argument('-d', '--data_directory_path', required=True,
                            help='Path to data directory')
        parser.add_argument('-p', '--chr_ploidy_file_path', required=True, type=str,
                            help="Path to chromosome ploidy file")
        parser.add_argument('-r', '--reference_genome_file_path', required=True,
                            help="Fasta file containing the reference genome")
        parser.add_argument('-i', '--ignore_indels', action='store_true',
                            help="Use the reference genome base in place of an indel. "
                                 "Defaults to false.")
        parser.add_argument('-j', '--ignore_snps', action='store_true',
                            help="Use the reference genome base in place of a snp. "
                                 "Defaults to False.")
        parser.add_argument('--sample', default=None,
                            help='String representation of a Sample object. Must provide '
                                 'this argument or the "--sample_id" and "--gender" arguments.')
        parser.add_argument('-s', '--sample_id', type=int, default=None,
                            help='sample name in vcf when prepended with sample. Overrides id from '
                                 'the "--sample" argument, if it is provided.')
        parser.add_argument('-x', '--gender', action='store', choices=['male', 'female'], default=None,
                            help='Gender of input sample (male or female). Overrides gender from '
                                 'the "--sample" argument, if it is provided.')
        parser.add_argument('-c', '--chromosomes', type=lambda chrs: [chr_ for chr_ in chrs.split(',')],
                            help="optional, comma-separated chromosome list")

        args = parser.parse_args()

        config_parameters = {"ignore_snps": args.ignore_snps,
                             "ignore_indels": args.ignore_indels}
        genome_builder = GenomeBuilderStep(args.log_directory_path,
                                           args.data_directory_path,
                                           config_parameters)

        if args.sample:
            sample = eval(args.sample)
        else:
            #Create dummy sample for debug purposes.
            sample = Sample(None, "debug sample", None, None, None)

        #Update sample if sample_id and gender, if specified.
        if args.sample_id:
            sample.sample_id = args.sample_id
        if args.gender:
            sample.gender = args.gender

        reference_genome = CampareeUtils.create_genome(args.reference_genome_file_path)
        chr_ploidy_data = CampareeUtils.create_chr_ploidy_data(args.chr_ploidy_file_path)
        genome_builder.execute(sample=sample,
                               chr_ploidy_data=chr_ploidy_data,
                               reference_genome=reference_genome,
                               chromosome_list=args.chromosomes)


class Genome:
    """
    Holds name, chromosome, current seq, current position (0 indexed) and current offset for a nascent, custom genome.
    The current offset is such that when it is added to the current position, one arrives at the corresponding position
    (0 indexed) on the reference genome.  The object also provides methods for appending, inserting and deleting based
    upon instructions in the variants input file.
    """

    #Patterns used to name all of the output files generated by the Genome class:

    #FASTA file containing the custom genome's sequence
    GENOME_OUTPUT_FILENAME_PATTERN = ('{genome_output_file_stem}_{genome_name}.fa')
    #List of all indel operations performed on the reference, to generate the
    #custom genome.
    INDEL_OUTPUT_FILENAME_PATTERN = ('{genome_output_file_stem}_indels_{genome_name}.txt')
    #File mapping coordinates from custom genome to their corresponding coordinates
    #in the reference.
    GENOME_MAPPER_FILENAME_PATTERN = ('{genome_output_file_stem}_to_parental_coordinate_mapping_{genome_name}.txt')

    def __init__(self, name, chromosome, start_sequence, start_position, genome_output_file_stem):
        self.name = name
        self.chromosome = chromosome
        self.sequence = StringIO()
        self.sequence.write(start_sequence)
        self.position = start_position
        self.offset = 0
        self.genome_output_filename = \
            self.GENOME_OUTPUT_FILENAME_PATTERN.format(
                genome_output_file_stem=genome_output_file_stem,
                genome_name=self.name
            )
        self.genome_indels_filename = \
            self.INDEL_OUTPUT_FILENAME_PATTERN.format(
                genome_output_file_stem=genome_output_file_stem,
                genome_name=self.name
            )
        self.indels_file = open(self.genome_indels_filename, 'a')
        self.genome_mapper_filename = \
            self.GENOME_MAPPER_FILENAME_PATTERN.format(
                genome_output_file_stem=genome_output_file_stem,
                genome_name=self.name
            )
        self.mapper_file = open(self.genome_mapper_filename, 'a')

        # These variables are used to control the mapper.  The last construct indicates what the last construct event
        # was (A = append, I = insert, D = delete).  The placeholders for the reference and this genome indicate where
        # the last append after an insert or delete started so we can mash all consecutive appends (due to snps)
        # together.
        #self.mapper_file.write("#CHR\tREF GENOME\tNEW GENOME\n")
        self.last_construct = "A"
        self.ref_placeholder = 0
        self.gen_placeholder = 0

    def append_segment(self, sequence):
        """
        Append the given sequence segment to the custom genome.  Since the sequence segment either has a one to one
        correspondence with that of reference genome or is a sequence segment drawn from the reference genome;
        execution of this method does not alter the current position of the custom genome relative to the current
        position of the reference genome/variant.  So position advances by the sequence segment length but offset
        remains unchangeed.
        :param sequence: sequence segment to append
        """

        # Last construct now was an append
        self.last_construct = "A"

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

        # Spans just before this insert (only done if the last construct was an append)
        if self.last_construct == "A":

            # Starts bumped for 1 index, ends not bumped because we are 1 base prior to insert
            reference_span = f"{self.ref_placeholder + 1}-{self.position + self.offset}"
            genome_span = f"{self.gen_placeholder + 1}-{self.position}"
            self.mapper_file.write(f"{self.chromosome}\t{reference_span}\t{genome_span}\n")

        # Span for this insert
        # Genome start bumped for 1 index but end not bumped because length includes starting base.
        genome_span = f"{self.position + 1}-{self.position + len(sequence)}"
        self.mapper_file.write(f"{self.chromosome}\t*-*\t{genome_span}\n")

        # Placeholders set to last positions after insert
        self.ref_placeholder = self.position + self.offset
        self.gen_placeholder = self.position + len(sequence)

        # Last construct now was an insert
        self.last_construct = "I"

        # This does the real insert work
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

        # Spans just before this delete (only done if the last construct was an append)
        if self.last_construct == "A":

            # Starts bumped for 1 index, ends not bumped because we are 1 base prior to delete
            reference_span = f"{self.ref_placeholder + 1}-{self.position + self.offset}"
            genome_span = f"{self.gen_placeholder + 1}-{self.position}"
            self.mapper_file.write(f"{self.chromosome}\t{reference_span}\t{genome_span}\n")

        # Span for this delete
        # Ref start bumped for 1 index but end not bumped because length includes starting base
        reference_span = f"{self.position + self.offset + 1}-{self.position + self.offset + length}"
        self.mapper_file.write(f"{self.chromosome}\t{reference_span}\t*-*\n")

        # Placeholders set to last positions after delete
        self.ref_placeholder = self.position + self.offset + length
        self.gen_placeholder = self.position

        # Last construct now was a delete
        self.last_construct = "D"

        # This does the real delete work
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
        self.mapper_file.close()

    def __str__(self):
        """
        Provide a string representation of the Genome object mainly for debugging purposes.
        :return: string representation of the Genome object
        """
        return f"name: {self.name}, chromosome: {self.chromosome}, position: {self.position}, offset: {self.offset}"


if __name__ == "__main__":
    sys.exit(GenomeBuilderStep.main())

    '''
    #The code below was originally here for diagnostic purposes. Command line input
    #is now processed directly by the main() method, and is used to submit the
    #GenomeBuilderStep to the system job scheduler (e.g. LSF, SGE).

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
-d ../../data/_run600/expression_pipeline/data \
-r ../resources/index_files/baby_genome.mm10/baby_genome.mm10.oneline_seqs.fa \
-p ../resources/index_files/baby_genome.mm10/baby_genome.mm10.chr_ploidy.txt \
-l ../../data/_run600/expression_pipeline/logs \
-x male
'''
