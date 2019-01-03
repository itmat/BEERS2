import re
from io import StringIO
import gzip

class ExpressionUtils:

    @staticmethod
    def edit_reference_genome(reference_genome_filename, edited_reference_genome_filename):
        reference_genome = dict()
        fasta_chromosome_pattern = re.compile(">([^\s]*)")
        chromosome, sequence = '', None
        building_sequence = False
        with open(reference_genome_filename, 'r') as reference_genome_file:
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
        with open(edited_reference_genome_filename, 'w') as edited_reference_genome_file:
            for chr, seq in reference_genome.items():
                edited_reference_genome_file.write(f">{chr}\n")
                edited_reference_genome_file.write(f"{seq}\n")

    @staticmethod
    def create_reference_genome(reference_genome_file_path):
        """
        Creates a reference genome dictionary from the reference genome file located at the provided path
        (assumed to be uncompressed).  The filename is assumed to contain the chr sequences without line breaks.
        :param reference_genome_file_path: path to reference genome
        :return: reference genome as a dictionary with the chromosomes/contigs as keys and the sequences as values.
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
        return reference_genome

if __name__ == "__main__":
    ExpressionUtils.edit_reference_genome("../../resources/index_files/hg38/Homo_sapiens.GRCh38.reference_genome.fa", "../../resources/index_files/hg38/Homo_sapiens.GRCh38.reference_genome_edited.fa")
