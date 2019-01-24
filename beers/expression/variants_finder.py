import os
import pysam
import re
import sys
from collections import namedtuple
from operator import attrgetter, itemgetter
import math
from io import StringIO
from beers.constants import CONSTANTS
from prettytable import PrettyTable


Read = namedtuple('Read', ['position', 'description'])
Read.__doc__ = """
A named tuple that possesses all the attributes of a variant
type:  match (M), deletion (D), insertion (I)
chromosome: chrN
position: position on ref genome
description: description of the variant (e.g., C, IAA, D5, etc.) 
"""


class VariantsFinderStep:
    """
    This class creates a text file listing variants for those locations in the reference genome having variants.
    The variants include snps and indels with the number of reads attributed to each variant.
    The relevant bam-formatted input file is expected to be indexed and sorted.

    This script outputs a file that gives the full breakdown at each
    location in the genome of the number of A's, C's, G's and T's as
    well as the number of each size of insertion and deletion.
    If it's an insertion the sequence of the insertion itself is given.
    So for example a line of output like the following means
    29 reads had a C in that location and three reads had an insertion of TTT.
    chr1:10128503 | C:29 | ITTT:3

    Note that only the top two variants are kept and of those the lesser variant's counts must meet
    certain user criteria (minimum threshold, read total count) to be considered a variant.  Single
    reads that match the corresponding base in the reference genome are not variants and as such are not kept.
    """

    DEFAULT_MIN_THRESHOLD = 0.03

    name = "Variants Finder Step"

    def __init__(self, log_directory_path, data_directory_path, parameters = dict()):
        self.data_directory_path = data_directory_path
        self.entropy_sort = parameters.get("sort_by_entropy", False)
        self.min_abundance_threshold = parameters.get('min_threshold', VariantsFinderStep.DEFAULT_MIN_THRESHOLD)
        self.clip_at_start_pattern = re.compile(r"(^\d+)[SH]")
        self.clip_at_end_pattern = re.compile(r"\d+[SH]$")
        self.variant_pattern = re.compile(r"(\d+)([NMID])")
        self.indel_pattern = re.compile(r"\|([^|]+)")
        self.log_directory_path = log_directory_path

    def validate(self):
        valid = True
        if not (0 < self.min_abundance_threshold < 1):
            print(f"The min_threshold, {self.min_abundance_threshold} must be between 0 and 1 exclusive.",
                  file=sys.stderr)
            valid = False
        return valid

    def remove_clips(self, cigar, sequence):
        """
        Remove soft and hard clips at the beginning and end of the cigar string and remove soft and hard clips at
        the beginning of the seq as well.  Modified cigar string and sequence are returned
        :param cigar: raw cigar string from read
        :param sequence: raw sequence string from read
        :return: tuple of modified cigar and sequence strings (sans clips)
        """
        clip_at_start = re.search(self.clip_at_start_pattern, cigar)
        if clip_at_start:
            cigar = re.sub(self.clip_at_start_pattern, "", cigar)
            sequence = sequence[int(clip_at_start.group(1)):]
        cigar = re.sub(self.clip_at_end_pattern, "", cigar)
        return cigar, sequence

    def call_variants(self, chromosome, reads):
        """
        Parses the reads dictionary (read named tuple:read count) for each chromosome - position to create
        a line with the variants and their counts delimited by pipes.  Dumping each chromosome's worth of
        data at a time is done to avoid too sizable a dictionary.  Additionally, if the user requests a sort by entropy,
        this function will do that ordering and send that data to stdout.
        :param chromosome: chromosome under consideration here
        :param reads: dictionary of reads to read counts
        """

        # variants list
        variants = []

        # Initializing the variable that will hold position information objects
        position_info = None

        # This dictionary is only used if the user requests that the read lines be sorted by entropy
        entropy_map = dict()

        # Iterate over the reads in the dictionary of variants to read counts sorted by the read position
        for read in sorted(reads.keys(), key=attrgetter('position')):

            # Initial iteration - set up position information object.
            if not position_info:
                position_info = PositionInfo(chromosome, read.position)

            # If the new position differs from the position of the position information currently being
            # consolidated, dump the current position information to the variants file if it is determined to
            # contain at least one variant.  In either case, create a new position information object for the new
            # position.
            if read.position != position_info.position:
                self.identify_variant(position_info, variants)

                # If the sort by entropy option is selected, also add to the entropy map dictionary the position
                # information entropy, keyed by the line content but only if the total number of reads exceeds the
                # depth cutoff.
                if self.entropy_sort and position_info.get_total_reads() >= int(self.depth_cutoff):
                        entropy_map[position_info.__str__()] = position_info.calculate_entropy()

                position_info = PositionInfo(chromosome, read.position)

            # Add the read description and read count to the position information
            position_info.add_read(read.description, reads[read])

        # Now that the reads are exhausted for this chromosome, dump the data from the current position information
        # object to the file.  Note that there may be no variants.
        self.identify_variant(position_info, variants)

        # If the user selected the sort by entropy option, other the entropy_map entries in descending order
        # of entropy and print to std out.
        if self.entropy_sort:
            sorted_entropies = sorted(entropy_map.items(), key=itemgetter(1), reverse=True)
            for key, value in sorted_entropies:
                print(key, end='')
        return variants

    def identify_variant(self, position_info, variants):
        """
        Helper method to filter position reads to identify variants
        :param position_info: position being evaluated
        :param variants: growing list of variants to which this position may be added if it contains variants
        """
        if position_info:
            reference_base = self.reference_genome[position_info.chromosome][position_info.position - 1]
            position_info.filter_reads(self.min_abundance_threshold, reference_base)
            if position_info:
                variants.append(position_info)

    def collect_reads(self, chromosome):
        """
        Iterate over the input txt file containing cigar, seq, start location, chromosome for each read and consolidate
        reads for each position on the genome.
        """

        reads = dict()
        current_read_components = []
        current_start = None

        for line in self.alignment_file.fetch(chromosome):

            # Remove unaligned reads, reverse reads, and non-unique alignments
            if line.is_unmapped or not line.is_read1 or line.get_tag(tag="NH") != 1:
                continue

            # Alignment Seqment reference_start is zero-based - so adding 1 to conform to convention.
            start = line.reference_start + 1
            sequence = line.query_sequence
            cigar = line.cigarstring
            cigar, sequence = self.remove_clips(cigar, sequence)
            current_pos_in_genome = int(start)
            loc_on_read = 1

            # Dropping duplicate reads (assumed to be PCR artifacts)
            if not current_start:
                current_start = start
                current_read_components.append((cigar, sequence))
            elif start != current_start:
                current_start = start
                current_read_components = []
                current_read_components.append((cigar, sequence))
            elif (cigar, sequence) in current_read_components:
                continue
            else:
                current_read_components.append((cigar, sequence))

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
                        key = Read(location, sequence[loc_on_read - 1])
                        reads[key] = reads.get(key, 0) + 1
                        loc_on_read += 1
                        current_pos_in_genome += 1
                    continue

                # For a deletion, designate the read named tuple description with a Dn where n is the
                # length of the deletion starting at this position.  In this way, subsequent reads having a
                # deletion of the same length at the same position will be added to this key.
                if read_type == "D":
                    location = current_pos_in_genome
                    key = Read(location, f'D{length}')
                    reads[key] = reads.get(key, 0) + 1
                    current_pos_in_genome += length
                    continue

                # For an insert, designate the read named tuple description with an Ib+ where b+ are the
                # bases to a inserted starting with this position.  In this way, subsequent reads having an
                # insertion of the same bases at the same position will be added to this key.
                if read_type == "I":
                    location = current_pos_in_genome
                    insertion_sequence = sequence[loc_on_read - 1: loc_on_read - 1 + length]
                    key = Read(location, f'I{insertion_sequence}')
                    reads[key] = reads.get(key, 0) + 1
                    loc_on_read += length
        return reads

    def execute(self, sample, alignment_file_path, chr_ploidy_data, reference_genome, chromosomes = None):
        """
        Entry point into variants_finder.  Iterates over the chromosomes in the list provided by the chr_ploidy_data
        keys to pick out variants.  Chromosomes that are not pertainent to the sample's gender are skipped.  If no
        sample gender is specified, only those chromosomes that have the same ploidy for both genders are processed.
        :param sample: The sample for which the variants for to be found
        :param chr_ploidy_data: dictionary of chromosomes as keys and a dictionary of male/female ploidy as values.
        :param reference_genome: A dictionary representation of the reference genome
        :param chromosomes: A listing of chromosomes to replace the list obtained from the alignment file.  Used for
        debugging purposes.
        """
        variants_filename = CONSTANTS.VARIANTS_FILE_NAME
        variants_file_path = os.path.join(self.data_directory_path, f'sample{sample.sample_id}', variants_filename)
        log_file_path = os.path.join(self.log_directory_path, f'sample{sample.sample_id}', __class__.__name__ + ".log")
        self.alignment_file = pysam.AlignmentFile(alignment_file_path, "rb")
        self.chromosomes = chromosomes if chromosomes else chr_ploidy_data.keys()
        self.reference_genome = reference_genome
        self.filter_chromosome_list(sample, chr_ploidy_data)
        log_table = PrettyTable()
        log_table.field_names =['chromosome','chromosome length','# positions with variants',
                                '# variants having no ref base variant','# positions having 1 variant',
                                '# positions having 2 variants']
        log_table.align['chromosome length'] = 'r'
        log_table.align['# positions with variants'] = 'r'
        log_table.align['# positions having no ref base variant'] = 'r'
        log_table.align['# positions having 1 variant'] = 'r'
        log_table.align['# positions having 2 variants'] = 'r'
        row_totals = [0, 0, 0, 0, 0]
        for chromosome in self.chromosomes:
            print(f"Finding variants for chromosome {chromosome}")
            variants = self.call_variants(chromosome, self.collect_reads(chromosome))
            self.load_variants(variants, variants_file_path)
            variants_without_ref_base = len([variant for variant in variants if not variant.contains_reference_base])
            pos_with_one_variant = len([variant for variant in variants if len(variant.reads) == 1])
            pos_with_two_variants = len([variant for variant in variants if len(variant.reads) == 2])
            row_values = [len(self.reference_genome[chromosome]), len(variants),
                          variants_without_ref_base, pos_with_one_variant, pos_with_two_variants]
            row_totals = [sum(item) for item in zip(row_totals, row_values)]
            row_values = [chromosome] + row_values
            log_table.add_row(row_values)
        row_totals = ['Totals'] + row_totals
        log_table.add_row(row_totals)
        with open(log_file_path, 'w') as log_file:
            log_file.write(log_table.get_string())
            log_file.write('\n')

    def filter_chromosome_list(self, sample, chr_ploidy_data):
        """
        Culls from the chromosome list, those chromosomes that are either not relevant given the sample gender or
        not relevant because no sample gender was provided.
        :param sample: subject sample which contains gender information
        :param chr_ploidy_data: dictionary of chromosomes as keys and a dictionary of male/female ploidy as values.
        """
        gender = sample.gender
        if not gender:
            self.chromosomes = [chr_ for chr_ in self.chromosomes
                                if chr_ploidy_data[chr_][CONSTANTS.MALE_GENDER] ==
                                chr_ploidy_data[chr_][CONSTANTS.FEMALE_GENDER] != 0]
        else:
            self.chromosomes = [chr_ for chr_ in self.chromosomes if chr_ploidy_data[chr_][gender] != 0]

    def load_variants(self, variants, variants_file_path):
        """
        Load the variants to a file in the user's designated output directory one chromosome at a time.  The filename
        has the stem of the alignment filename suffixed with _variants.txt
        :param variants: variants list for one chromosome.
        """
        with open(variants_file_path, 'a') as variants_file:
            for variant in variants:
                variants_file.write(variant.__str__())


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
        self.contains_reference_base = None

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

    def filter_reads(self, min_abundance_threshold, reference_base):
        """
        Filters out from this position, reads that are not considered true variants.  Any reads with read counts of
        only 1 are excluded to start with. At most, only the top two remaining reads are retained.  The lesser of those
        two reads may also be removed if it does not satisfy the minimum abundance threshold criterion.  The minimum
        abundance threshold criterion specifies that the percent contribution of the lesser variant reads to the total
        reads be equal or greater than the threshold provided.  In the event of a tie for one of both of those top
        two slots, preference is given to the reference base if it is included in the tie. If at any point in
        filtering, only one read remains and its description matches the reference base, it is removed, leaving no
        variants.  Once complete, the reads for this position object contain only true variants (which may include
        the reference base if there is one another true variant).
        :param min_abundance_threshold:  criterion for minimum abundance threshold
        :param reference_base:  the base of the reference genome at this position.
        """

        variants = []

        # Remove all reads with read counts of no more than 1
        filtered_reads = [read for read in self.reads if read[1] > 1]

        # If only a single read remains, remove if it matches the reference base - there are no variants
        if len(filtered_reads) == 1 and filtered_reads[0][0] == reference_base:
            variants = []

        # If multiple reads remain, retain only the top 2 reads
        elif len(filtered_reads) > 1:

            candidate_variants = {filtered_read[0]: filtered_read[1] for filtered_read in filtered_reads}

            for _ in range(2):
                max_read = max(candidate_variants.items(), key=itemgetter(1))[0]
                variants.append((max_read, candidate_variants[max_read]))
                del candidate_variants[max_read]

            # Get reference read, if it exists
            reference_read = [(description, read_count) for description, read_count in filtered_reads
                              if description == reference_base]

            # If the reference read exists and is not among the variants but has the same number of counts as
            # the second place read, replace the second place read with it.
            if reference_read and reference_read[0] not in variants and reference_read[0][1] == variants[-1][1]:
                variants[-1] = reference_read[0]

            # Remove the second place read if it does not satisfy the min_abundance_threshold criterion
            total_reads = sum(read_count for _, read_count in variants)
            if variants[1][1] / total_reads < min_abundance_threshold:
                if variants[0][0] != reference_base:
                    variants = [variants[0]]
                else:
                    variants = []

        # The only position reads remaining are variants
        if variants:
            self.contains_reference_base = any([variant[0] == reference_base for variant in variants])
        self.reads = variants

    def __bool__(self):
        """
        The object's truth value depends on the presence of reads.  Following a filter_reads step, that value will
        essentially depend on the presence of variants.
        :return: True if reads (variants) are present and false otherwise.
        """
        return True if self.reads else False

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
