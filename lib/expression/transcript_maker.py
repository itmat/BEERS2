import sys
import argparse
import re
import os
import numpy as np
import bisect
sys.path.append('/Users/crislawrence/Documents/Work/BEERS/BEERS2.0/lib/beers')
from molecule import Molecule
# from sample import Sample


class TranscriptMaker:

    def __init__(self, annotation_filename,
                 quantities_filename,
                 genome_filename,
                 log_filename,
                 seed,
                 transcriptome_size):
        self.quantities_filename = quantities_filename
        self.annotation_filename = annotation_filename
        self.genome_filename = genome_filename
        self.log_filename = log_filename

        # TODO - create and log a seed if none provided.
        self.seed = int(seed)
        np.random.seed(self.seed)

        self.transcriptome_size = int(transcriptome_size)
        self.polya_tail = 200 * 'A'

        # Holds molecules created from transcript
        self.molecules = set()

        # Holds the cumulative transcript distribution
        self.transcript_distribution = TranscriptCumulativeDistribution()

        # Hold the probabalistically determined transcript counts
        self.transcript_counts_map = dict()

        # Holds unique listing of exon locations
        self.exon_location_list = set()

        # Provides the regex pattern for extracting components of the exon location string
        self.exon_info_pattern = re.compile('(.*):(\d+)-(\d+)')

        # Dictionaries to record whether a chromosome is available in the genome file, is available in the exon
        # file.
        self.chromosome_in_genome_file = dict()
        self.chromosome_in_exon_file = dict()

        # Since the log file will be open for appending, we need to insure that the file doesn't
        # currently exist.
        try:
            os.remove(self.log_filename)
        except OSError:
            pass

    def create_transcript_cumulative_distribution(self):
        with open(self.quantities_filename, 'r') as quantities_file:
            for line in quantities_file:
                if line.startswith(' '):
                    continue

                fields = line.rstrip('\n').split('\t')

                # Transcript id and counts
                transcript_id, transcript_counts = fields[0], int(fields[1])
                self.transcript_distribution.append(transcript_id, transcript_counts)

    def create_transcript_counts_map(self):
        for _ in range(self.transcriptome_size):
            transcript_id = self.transcript_distribution.locate_id(np.random.uniform())
            self.transcript_counts_map[transcript_id] = self.transcript_counts_map.get(transcript_id, 0) + 1

    def prepare_transcripts(self):

        # Create a map of the transcript distributions
        self.create_transcript_cumulative_distribution()
        self.transcript_distribution.normalize()

        # Pair transcript_id with counts probabilistically
        self.create_transcript_counts_map()
        print(sum(self.transcript_counts_map.values()))

        # Create a unique listing of exon locations from the transcript annotation file data.
        self.create_exon_location_list()

        # Open the genome fasta file for reading only.
        with open(self.genome_filename, 'r') as genome_file:

            # Iterate over each chromosome in the genome fasta file.  Note that the chromosome sequence is expected on
            # only one line at this stage.
            for line in genome_file:

                # Remove the leading '>' character to leave the chromosome
                chromosome = line.lstrip('>').rstrip('\n')

                # TODO - awkward hack to compare UCSC and ensembl
                chromosome = chromosome[3:]
                if chromosome == 'M':
                    chromosome += 'T'

                # Note that the chromosome 'chromosome' is among those listed in the genome file
                self.chromosome_in_genome_file[chromosome] = True

                # Collect the chromosome sequence from the following line
                sequence = genome_file.readline().rstrip('\n')

                # Use the chromosome and its sequence to construct a mapping of exon location string to
                # their sequences.
                exon_sequence_map = self.create_exon_sequence_map(chromosome, sequence)
                print(f"Done with exons for {chromosome}")

                # Generate the gene fasta file using the genome chromosome and the related exon sequence map.
                total_per_chrom = self.make_molecules(chromosome, exon_sequence_map)
                print(f"Done with {total_per_chrom} molecules for {chromosome}")

    def create_exon_location_list(self):
        """
        Generate a unique listing of exon location strings from the provided transcript annotation file.  Note that
        the same exon may appear in multiple transcripts.  So the listing is actually a set to avoid duplicate entries.
        """

        # Open the transcript annotation file for reading only
        with open(self.annotation_filename, 'r') as annotation_file:

            # Iterate over each line in the file
            for line in annotation_file:

                # Collect all the field values for the line read following newline removal.
                chromosome, strand, start, end, exon_count, exon_starts, exon_ends, transcript_id, *other = \
                    line.rstrip('\n').split('\t')

                # TODO Ick!  Need a better way to tell if there is a header
                if chromosome == 'chrom':
                    continue

                # Remove any trailing commas in the exon starts and exon ends fields and split
                # the starts and stops into their corresponding lists
                exon_starts_list = re.sub(r'\s*,\s*$', '', exon_starts).split(",")
                exon_ends_list = re.sub(r'\s*,\s*$', '', exon_ends).split(",")

                # For each exon belonging to the transcript, construct the exon's location string and add it to the
                # growing set of exon locations assuming it is not already present.
                for index in range(int(exon_count)):
                    exon_location = f'{chromosome}:{(int(exon_starts_list[index]) + 1)}-{(exon_ends_list[index])}'
                    self.exon_location_list.add(exon_location)

    def create_exon_sequence_map(self, genome_chromosome, sequence):
        """
        For the given genome chromosome and its sequence, create a dictionary of exon sequences keyed to the exon's
        location (i.e., chr:start-end).
        :param genome_chromosome: given genome chromosome
        :param sequence: the genome sequence corresponding to the genome chromosome (without line breaks)
        :return: map of exon location : exon sequence
        """

        # Start with a empty dictionary
        exon_sequence_map = dict()

        # Iterate over the list of previously obtained exon locations
        for exon_location in self.exon_location_list:

            # Extract the chromosome, start and end from the exon location string
            exon_info_match = re.search(self.exon_info_pattern, exon_location)
            chromosome = exon_info_match.group(1)
            exon_start = int(exon_info_match.group(2))
            exon_end = int(exon_info_match.group(3))

            # Note that the exon listing contains at least one exon on chromosome 'chromosome'
            self.chromosome_in_exon_file[chromosome] = True

            # If the chromosome on which the exon is located is the same of the genome chromosome
            # provided as a parameter, get the sequence for that exon and create a dictionary
            # entry relating the exon location string to the exon sequence.
            if chromosome == genome_chromosome:
                exon_sequence_map[exon_location] = sequence[exon_start - 1:exon_end]

        # Return the dictionary of exon location string : exon sequence for the genome chromosome
        # provided in the parameter list.
        return exon_sequence_map

    def make_molecules(self, genome_chromosome, exon_sequence_map):

        total = 0
        # Open the transcript annotation file for reading
        with open(self.annotation_filename, 'r') as annotation_file, open(self.log_filename, 'a') as log_file:

                # Iterate over the annotation file
                for line in annotation_file:

                    # Collect all the field values for the line read following newline removal.
                    chromosome, strand, start, end, exon_count, exon_starts, exon_ends, transcript_id, *other = \
                        line.rstrip('\n').split('\t')

                    # TODO Ick!  Need a better way to tell if there is a header
                    if chromosome == 'chrom':
                        continue

                    # Get quantity for this transcript
                    quantity = self.transcript_counts_map.get(transcript_id, 0)

                    if quantity > 0:

                        # Remove trailing commas in the exon starts and exon ends fields and split
                        # the starts and stops into their corresponding lists
                        exon_starts_list = re.sub(r'\s*,\s*$', '', exon_starts).split(",")
                        exon_ends_list = re.sub(r'\s*,\s*$', '', exon_ends).split(",")

                        # If the chromosome on which this transcript is located is the same as the genome chromosome
                        # provided as a parameter, render this transcript as molecules, the number given by the
                        # transcript counts.
                        if chromosome == genome_chromosome:

                            # Initialize the gene's sequence
                            transcript_sequence = ""

                            # For each exon belonging to the transcript, construct the exon's location string and use it as
                            # a key to obtain the actual exon sequence.  Concatenate that exon sequence to the transcript
                            # sequence.  Note that the 1 added to each exon start takes into account the zero based
                            # and half-open ucsc coordinates.
                            for index in range(int(exon_count)):
                                exon_key = f'{chromosome}:{int(exon_starts_list[index]) + 1}-{(exon_ends_list[index])}'
                                exon_sequence = exon_sequence_map[exon_key]
                                transcript_sequence += exon_sequence

                            # Create an initial cigar sequence (all matches)
                            cigar = str(len(transcript_sequence)) + 'M'

                            # Log the transcript
                            log_file.write(f'{transcript_id}\t{quantity}\t{transcript_sequence}\n')
                            total += quantity

                            # Affix a poly A tail
                            transcript_sequence += self.polya_tail

                            # Generate as many molecules of the transcript as dictated by the transcript's quantity.
                            for _ in range(quantity):
                                self.molecules.add(Molecule(Molecule.next_molecule_id, transcript_sequence, 1,
                                                            cigar, transcript_id))
                                Molecule.next_molecule_id += 1
        return(total)

    @staticmethod
    def main():
        parser = argparse.ArgumentParser(description='Transcript_Maker')
        parser.add_argument('-a', '--maternal_annot_filename')
        parser.add_argument('-q', '--maternal_quant_filename')
        parser.add_argument('-g', '--maternal_genome_filename')
        parser.add_argument('-x', '--paternal_annot_filename')
        parser.add_argument('-y', '--paternal_quant_filename')
        parser.add_argument('-z', '--paternal_genome_filename')
        parser.add_argument('-l', '--log_filename')
        parser.add_argument('-d', '--seed')
        parser.add_argument('-s', '--transcriptome_size')
        args = parser.parse_args()
        print(args)
        sample_molecules = set()
        # transcript_maker1 = TranscriptMaker(args.maternal_annot_filename,
        #                                    args.maternal_quant_filename,
        #                                    args.maternal_genome_filename,
        #                                    args.log_filename,
        #                                    args.transcriptome_size)
        # transcript_maker1.prepare_transcripts()
        # sample_molecules.union(transcript_maker1.molecules)
        transcript_maker2 = TranscriptMaker(args.paternal_annot_filename,
                                            args.paternal_quant_filename,
                                            args.paternal_genome_filename,
                                            args.log_filename,
                                            args.seed,
                                            args.transcriptome_size)
        transcript_maker2.prepare_transcripts()
        sample_molecules.union(transcript_maker2.molecules)
        print(f'Number of molecules: {len(transcript_maker2.molecules)}')
        # return Sample(sample_molecules)


class TranscriptCumulativeDistribution:

    def __init__(self):
        self.transcript_ids = []
        self.transcript_cumulative_counts = []
        self.count_accumulator = 0

    def append(self, transcript_id, count):
        self.transcript_ids.append(transcript_id)
        self.count_accumulator += count
        self.transcript_cumulative_counts.append(self.count_accumulator)

    def normalize(self):
        self.transcript_cumulative_counts = \
            [count/self.count_accumulator for count in self.transcript_cumulative_counts]

    def locate_id(self, value):
        return self.transcript_ids[bisect.bisect(self.transcript_cumulative_counts, value)]


if __name__ == "__main__":
    sys.exit(TranscriptMaker.main())

'''
Example command
python transcript_maker.py \
    -a ../../data/preBEERS/mm10/Mus_musculus.GRCm38.93.annotation.txt \
    -q ../../data/preBEERS/mm10/m_tx_quant.tsv \
    -g ../../data/preBEERS/mm10/reference_genome_edited.fa \
    -x ../../data/preBEERS/mm10/Mus_musculus.GRCm38.93.annotation.txt \
    -y ../../data/preBEERS/mm10/p_tx_quant.tsv \
    -z ../../data/preBEERS/mm10/reference_genome_edited.fa \
    -l ../../data/preBEERS/mm10/transcript_maker.log \
    -d 100 \
    -s 1000000
    '''
