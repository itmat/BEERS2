from io import StringIO
import numpy as np
import sys
import re
import os
import argparse
from timeit import default_timer as timer
from operator import itemgetter
import time


class GenomeMaker:
    """
    Used to create two custom genomes based upon the contents of a variants file and an appropriate reference genome
    file.  The custom genomes are built one chromosome at a time.
    """

    def __init__(self,
                 variants_filename,
                 genome_ref_filename,
                 genome_output_file_stem,
                 seed,
                 min_threshold,
                 log_filename,
                 gender,
                 gender_chr_names,
                 min_read_total_count,
                 ignore_indels,
                 ignore_snps):
        """
        Constructor for the GenomeMaker object which holds the arguments provided.  It also contains the hardcoded
        custom genome names and the regular expression used to extract the needed data from the variants file.  Since
        the output files holding the custom genomes are constructed one chromosome at a time, they must start out as
        virgin files.  So any prior files having the same path/filename are first deleted.
        :param variants_filename: path of the text file containing the variants (formatted as in the variants_maker.py
        file in this package).  Below is an example of the expected variants file formatting:
            chr1:1498605 | T:12 | D1:1      TOT=13  0.92,0.08       E=0.39124356362925566
            chr1:1498606 | T:13 | IT:1      TOT=14  0.93,0.07       E=0.37123232664087563
        :param genome_ref_filename:  fasta file containing the reference genome (each chromosome sequence should be on
        a single line - no line breaks)  Also note that the chromosome order must agree with the order of the
        chromosomes in the variants_filename because once passed, as reference chromosome is not revisited.
        :param genome_output_file_stem: path to the output fasta file(s) - will be suffixed with the custom genome
        names.
        :param seed: seed for random number generator to allow for repeatable analyses.  Defaults to None
        :param min_threshold: minimum percent abundance of a variant for it to be recognized as a legitimate variant.
        Defaults to 0.03
        :param log_filename: path to log file.
        :param gender: gender ascribed to the variants input file.
        :param gender_chr_names: listing of X, Y, M chromosome names in this order.
        :param min_read_total_count: when the total reads of the two most abundant variants meet or exceed this value,
        the second most abundant variant read count must exceed one to be included as a viable variant irregardless of
        whether it meets the minimum threshold requirement.
        :param ignore_indels: choose reference genome base over indel option when a indel variant is encountered.
        :param ignore_snps: choose reference genome base over snp option when a snp variant is encountered.
        """
        self.variants_filename = variants_filename
        self.genome_output_file_stem = genome_output_file_stem
        self.genome_ref_filename = genome_ref_filename
        self.log_filename = log_filename
        self.gender = gender
        self.gender_chr_names = {name[0]:name[1] for name in zip(['X','Y','M'],gender_chr_names)}
        self.ignore_indels = ignore_indels
        self.ignore_snps = ignore_snps

        # If seed is set, use it to assure reproducible results.  Otherwise generate a seed from a timestamp so
        # that it may be added to the log file in the event that the user wishes to repeat an experiment after the
        # fact.
        if not seed:
            seed = GenomeMaker.generate_seed()
        self.seed = seed
        np.random.seed(seed)
        self.min_abundance_threshold = min_threshold
        self.min_read_total_count = min_read_total_count

        # Pattern used to extract chromosome, position and variants from each variants line.
        # The [^\] is there to insure that we avoid greedy matches to subsequent colons since colons
        # delimit chromosome and position and also variant description and read count, but a pipe
        # separates the chromosome:position from the variant description:read count portion.
        self.variant_line_pattern = re.compile('^([^|]+):(\d+) \| (.*)\tTOT')
        self.genome_names = ['maternal', 'paternal']

        # Make sure the filenames for the genomes are pristine.
        for genome_name in self.genome_names:
            try:
                os.remove(self.genome_output_file_stem + "_" + genome_name + ".fa")
            except OSError:
                pass

        # Make sure the filenames for the genome indel lists are pristine.
        for genome_name in self.genome_names:
            try:
                os.remove(self.genome_output_file_stem + "_indels_" + genome_name + ".txt")
            except OSError:
                pass

    @staticmethod
    def generate_seed():
        """
        Provides an integer of 32 bits or less using a seconds based timestamp.  If the timestamp exceeds 32 bits, the
        integer obtained from the lowest 32 bits is returned
        :return: a 32 bit integer seed for the numpy random number generator.
        """
        candidate_seed = int(time.time())

        # If the timestamp bit length exceeds 32 bits, mask out all but the lowest 32 bits.
        if candidate_seed.bit_length() > 32:
            candidate_seed = candidate_seed & 0xffffffff
        return candidate_seed

    def get_most_abundant_variants(self, chromosome, variants):
        """
        Obtains the (at most) two variants for a chromosome and position having the most reads (greatest abundances).
        In many cases, there may be only one variant given, in which case it alone is returned.  In some cases, the
        lesser of the two most abundant variants may be so few reads as to be deemed invalid, in which case only the
        most abundant variant is returned.  Otherwise the top two variants are returned.
        :param chromosome: chromosome to which the variantes apply.  In some cases, only one chromosome of the given
        type is contributed, in which case, only the most abundant variant should be returned.
        :param variants: list of variant descriptions and their read counts
        :return: list of one or two more abundant variants sans read counts
        """
        max_variants = []

        # Make a fresh copy of the variants (since variants are mutable) that may be modified.  Create a dictionary
        # from the array (variant description: read count) and insure that reads are treated as integers.
        # (e.g., ['C:1','T:2','D2:2'] => {'C': 1, 'T': 2, 'D2': 2})
        variant_reads = {variant.split(':')[0]: int(variant.split(':')[1]) for variant in variants}

        # If only one variant is given, return it only without the read count.
        if len(variant_reads) == 1:
            return [next(iter(variant_reads.keys()))]

        # Otherwise, use the max operation to pluck off the top two variants, one at a time.
        for _ in range(2):

            # This line grabs the variant description with the maximum read count
            max_variant = max(variant_reads.items(), key=itemgetter(1))[0]

            # This line appends a tuple, (variant description, variant read count) to the max_variants arrays.
            max_variants.append((max_variant, variant_reads[max_variant]))

            # Remove the the variant_reads dictionary the entry having the max_variant as the key.
            del variant_reads[max_variant]

        # If the chromosome is unique to one parent (e.g., chrM or chrY) or unique due to gender (e.g. chr X and
        # gender is male) return the most abundant variant only.
        if chromosome in [self.gender_chr_names['M'], self.gender_chr_names['Y']] or (self.gender == 'male' and chromosome == self.gender_chr_names['X']):
            return [max_variants[0][0]]

        # Determine the total reads for the top two variants and if the lesser variant's percentage of the total
        # reads is below the threshold, discard it and return only the top variant without its read count.
        total_reads = sum(read for _, read in max_variants)
        if max_variants[1][1]/total_reads < self.min_abundance_threshold:
            return [max_variants[0][0]]

        # If the second most abundant variant has only 1 read when the total reads between the top two most abundant
        # variants is equal to or greater then the min_read_total_count, again return only the top variant without its
        # read count.
        if total_reads >= self.min_read_total_count and max_variants[1][1] == 1:
            return [max_variants[0][0]]

        # Otherwise, return both variants without the corresponding read counts.
        return [variant for variant, read in max_variants]

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

    def make_genome(self):
        """
        Performs the task of making two custom genomes, one chromosome at a time, from the provided variant file and
        genome reference file.
        """

        # Open the log file for writing
        with open(self.log_filename, 'w') as log_file:

            # Open both the variants file and genome reference file for reading
            with open(self.variants_filename) as variants_file, open(self.genome_ref_filename) as genome_ref_file:
                reference_chromosome = None
                reference_sequence = None

                # Flag to indicate when a genome (via the current chromosome) is being constructed.
                building_chromosome = False

                genomes = list()

                # Iterate over each variant line in the variants file
                for variant_line in variants_file:

                    # Extract the variant chromomsome, variant position and an array of the variants from the line.
                    match = re.match(self.variant_line_pattern, variant_line)
                    chromosome = match.group(1)

                    # Substraction to account for fact that variant text file indexed at 1
                    variant_position = int(match.group(2)) - 1
                    variants = match.group(3).split(' | ')

                    # Starting new chromosome in the variants file or, if no reference_chromosome yet exists, starting
                    # the very first chromosome.  We iterate in case the variants file skips over a chromosome.  The
                    # variants file and the reference genome fasta file must order the chromosomes in the same way but
                    # the reference genome file can be more complete (have more chromosomes) than the variants file.
                    while not reference_chromosome or reference_chromosome != chromosome:

                        # At least one chromosome already completed - finish and save
                        if building_chromosome:
                            for genome in genomes:
                                log_file.write(f"Appending"
                                               f" {len(reference_sequence[genome.position + genome.offset:])}"
                                               f" bases of reference sequence at reference position "
                                               f" {genome.position + genome.offset} to complete the genome for"
                                               f" chromosome {reference_chromosome}.\n")
                                genome.append_segment(reference_sequence[genome.position + genome.offset:])
                                log_file.write(f"{genome}\n")
                                genome.save_to_file()

                        # Get a reference sequence for the next chromosome in the reference genome file
                        line = genome_ref_file.readline()

                        # Accounts for the possibility that the reference genome file is exhausted before the
                        # variants file.  This can happen if the reference genome chromosomes and the variant
                        # chromosomes are not ordered identically.  This is a problem since both inputs are being
                        # streamed from files presently.
                        if not line:
                            print("Warning:  Variants file not exhausted but reference genome file is."
                                  "  Could be due to misordering.")
                            sys.exit(1)
                        if line.startswith(">"):
                            reference_chromosome = line[1:].rstrip()
                            print("Ref chr: " + reference_chromosome)

                            # Only set up genomes if the reference chromosome is the same at the one currently
                            # being read from the variants file.  Otherwise, skip over and get the next chromosome from
                            # the reference genome file.  An new set of of genomes is created for each new chromosome.
                            if reference_chromosome == chromosome:
                                building_chromosome = True
                                reference_sequence = genome_ref_file.readline().rstrip()
                                genomes = list()
                                log_file.write(f"Appending"
                                               f" {len(reference_sequence[:variant_position])}"
                                               f" bases of reference sequence at reference position "
                                               f" 0 to start both genomes for chromosome {chromosome}.\n")
                                start_sequence = reference_sequence[0: variant_position]

                                # In addition to diploid chromosomes, this genome will contain contributions that
                                # derive solely from the mother.  So this genome is skipped when building the Y
                                # chromosome.
                                if chromosome != self.gender_chr_names['Y']:
                                    genomes.append(
                                        Genome(self.genome_names[0], chromosome, start_sequence,
                                               variant_position, self.genome_output_file_stem))

                                # In addition to diploid chromosomes, this genome will contain contributions that
                                # derive solely from the father.  So this genome is skipped when building the M
                                # chromosome and when building the X chromosome if the input gender is male.  Also, it
                                # is possible that there may be Y variants even though the input gender is female.  If
                                # that is the case, neither genome is built for chromosome Y.
                                if chromosome != self.gender_chr_names['M'] \
                                        and not(self.gender == 'male' and chromosome == self.gender_chr_names['X'])\
                                        and not(self.gender == 'female' and chromosome == self.gender_chr_names['Y']):
                                    genomes.append(
                                        Genome(self.genome_names[1], chromosome, start_sequence,
                                               variant_position, self.genome_output_file_stem))

                                # There are no genomes to build move on to the next chromosome, if any.
                                if not genomes:
                                    building_chromosome = False
                                    continue
                            else:
                                genome_ref_file.readline()
                                building_chromosome = False

                    # Return the top two (based on number of reads) of the variants on this line.
                    max_variants = self.get_most_abundant_variants(chromosome, variants)

                    variant_keys = [variant.split(":")[0] for variant in variants]
                    if reference_sequence[variant_position] not in variant_keys:
                        log_file.write(f"Reference at {chromosome}:{variant_position} is "
                                       f"{reference_sequence[variant_position]} whereas variants "
                                       f"are {variant_keys}\n")

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

                            reference_base = reference_sequence[genome.position + genome.offset: variant_position]

                            # If only one variant, apply directly to genome
                            if len(max_variants) == 1:
                                self.build_sequence_from_variant(genome, max_variants[0], reference_base)

                            # If two variants exist, toss a coin for one
                            else:
                                max_variant = np.random.choice([max_variants[0], max_variants[1]], p=[0.5, 0.5])
                                self.build_sequence_from_variant(genome, max_variant, reference_base)
                                max_variants.remove(max_variant)


                # Variants file exhausted, must be done with last chromosome.  So save the last pair of genomes.
                for genome in genomes:
                    log_file.write(f"Appending"
                                   f" {len(reference_sequence[genome.position + genome.offset:])}"
                                   f" bases of reference sequence at reference position "
                                   f" {genome.position + genome.offset} to complete the genome for"
                                   f" chromosome {chromosome}.\n")
                    genome.append_segment(reference_sequence[genome.position + genome.offset:])
                    log_file.write(f"Final Genome for chromosome {reference_chromosome}: {genome}\n")
                    genome.save_to_file()

    @staticmethod
    def main():
        """
        Entry point for the genome maker.  Three input file paths/names are required.  One is the variant file, created
        by the variants maker in this package.  The second is a fasta file of the reference genome containing at least
        those chromosomes found in the variants file.  Each of the chromosome sequences in the fasta file must be
        resident on one line.  THe third is a log file.  One output file path/name is required, to which the hardcoded
        custom genome names will be suffixed to form fasta files of the custom genome sequences created.  Two optional
        arguments are the seed and the threshold.

        It should be noted that both the variants file and the genome reference fasta file should list the chromosomes
        in the same order.  Additionally, the variants file should list the variants in ascending order of position
        within each chromosome.
        """
        parser = argparse.ArgumentParser(description='Make Genome Files')
        parser.add_argument('-v', '--variants_filename',
                            help="Textfile providing variants as provided by the ouput"
                                 " of the variants_maker.py script.")
        parser.add_argument('-r', '--reference_genome_filename',
                            help="Fasta file containing the reference genome")
        parser.add_argument('-g', '--genome_output_file_stem',
                            help="Path of output genome file.  Will be suffixed _seq1.fa or _seq2.fa")
        parser.add_argument('-s', '--seed', type=int, default=None,
                            help="Integer to be used as a seed for the random number generator."
                                 "  Value defaults to no seed.")
        parser.add_argument('-m', '--min_threshold', type=float, default=0.03,
                            help="For any position, if the percent abundance of the second most abundant variant"
                                 "relative to the two most abundant variants falls below this threshold, the seond"
                                 "most abundant variant will be discarded.  Defaults to 0.03 (3 percent)")
        parser.add_argument('-l', '--log_filename', help="Log file.")
        parser.add_argument('-x', '--gender', action='store', choices=['male', 'female'],
                            help="Gender of input sample (male or female).")
        # TODO should be required - this will do under we have a pipeline
        parser.add_argument('-n', '--gender_chr_names', type=lambda names: [name for name in names.split(',')],
                            help="Enter the gender specific chromosome names in X, Y, M order.")
        parser.add_argument('-c', '--min_read_total_count', type=int, default=10,
                            help="For any position, if the total reads from the two most abundant variants"
                                 " equals or exceeds this value, the second most abundant variant will be discarded"
                                 " if it is a single read even if it passes the minimum threshold test.  Defaults to"
                                 " 10.")
        parser.add_argument('-i', '--ignore_indels', action='store_true',
                            help="Use the reference genome base in place of an indel.  Defaults to false.")
        parser.add_argument('-j', '--ignore_snps', action='store_true',
                            help="Use the reference genome base in place of a snp.  Defaults to False.")
        args = parser.parse_args()
        print(args)
        genome_maker = GenomeMaker(args.variants_filename,
                                   args.reference_genome_filename,
                                   args.genome_output_file_stem,
                                   args.seed,
                                   args.min_threshold,
                                   args.log_filename,
                                   args.gender,
                                   args.gender_chr_names,
                                   args.min_read_total_count,
                                   args.ignore_indels,
                                   args.ignore_snps)
        start = timer()
        genome_maker.make_genome()
        end = timer()
        sys.stderr.write(f"Genome Maker: {end - start} sec\n")


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
    sys.exit(GenomeMaker.main())


'''
Example Call:

python genome_maker.py \
-v  ../../data/expression/GRCh38/Test_data.1002_baseline.All_unique_mappers.fw_only.sorted_to_match_fasta_variants.txt \
-g ../../data/expression/GRCh38/genome_maker_output \
-r ../../data/expression/GRCh38/Homo_sapiens.GRCh38.reference_genome_edited.fa \
-l ../../data/expression/GRCh38/genome_maker.log -n "X,Y,MT" -s 100 -x female

'''
