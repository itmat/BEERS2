import re
from io import StringIO
import os
import gzip
import itertools

class ExpressionUtils:

    @staticmethod
    def edit_reference_genome(reference_genome_file_path, edited_reference_genome_file_path):
        """
        Helper method to convert a reference genome file containing line breaks embedded within its
        sequences to a reference genome file containing each seqeuence on a single line.
        :param reference_genome_file_path: Path to reference geneome file having multi-line sequence data
        :param edited_reference_genome_file_path: Path to reference genome file to create with single line sequence
        data.
        """
        reference_genome = dict()
        fasta_chromosome_pattern = re.compile(">([^\s]*)")
        chromosome, sequence = '', None
        building_sequence = False
        with open(reference_genome_file_path, 'r') as reference_genome_file:
            for line in reference_genome_file:
                if line.startswith(">"):
                    if building_sequence:
                        reference_genome[chromosome] = sequence.getvalue()
                        sequence.close()
                    chromosome_match = re.match(fasta_chromosome_pattern, line)
                    chromosome = chromosome_match.group(1)
                    building_sequence = True
                    sequence = StringIO()
                    continue
                elif building_sequence:
                    sequence.write(line.rstrip('\n').upper())
        with open(edited_reference_genome_file_path, 'w') as edited_reference_genome_file:
            for chr, seq in reference_genome.items():
                edited_reference_genome_file.write(f">{chr}\n")
                edited_reference_genome_file.write(f"{seq}\n")

    @staticmethod
    def create_reference_genome(reference_genome_file_path):
        """
        Creates a reference genome dictionary from the reference genome file located at the provided path
        (if compressed, it must have a gz extension).  The filename is assumed to contain the chr sequences without
        line breaks.
        :param reference_genome_file_path: path to reference genome file (either compressed or not)
        :return: reference genome as a dictionary with the chromosomes/contigs as keys and the sequences as values.
        """
        reference_genome = dict()
        fasta_chromosome_pattern = re.compile(">([^\s]*)")
        chromosome, sequence = '', None
        building_sequence = False
        _, file_extension = os.path.splitext(reference_genome_file_path)
        if 'gz' in file_extension:
            with gzip.open(reference_genome_file_path, 'r') as reference_genome_file:
                for chr, seq in itertools.zip_longest(*[reference_genome_file] * 2):
                    chr = chr.decode("ascii").rstrip()[1:]
                    reference_genome[chr] = seq.decode("ascii").rstrip()
        else:
            with open(reference_genome_file_path, 'r') as reference_genome_file:
                with open(reference_genome_file_path, 'r') as reference_genome_file:
                    for chr, seq in itertools.zip_longest(*[reference_genome_file] * 2):
                        chr = chr.rstrip()[1:]
                        reference_genome[chr] = seq.rstrip()
        return reference_genome

if __name__ == "__main__":
    ExpressionUtils.create_reference_genome("../../resources/index_files/hg38/Homo_sapiens.GRCh38.reference_genome.fa")
