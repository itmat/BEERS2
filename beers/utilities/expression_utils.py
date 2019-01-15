import re
from io import StringIO
import os
import gzip
import itertools
import pandas as pd

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
    def create_genome(genome_file_path):
        """
        Creates a genome dictionary from the genome file located at the provided path
        (if compressed, it must have a gz extension).  The filename is assumed to contain the chr sequences without
        line breaks.
        :param genome_file_path: path to reference genome file (either compressed or not)
        :return: genome as a dictionary with the chromosomes/contigs as keys and the sequences as values.
        """
        genome = dict()
        _, file_extension = os.path.splitext(genome_file_path)
        if 'gz' in file_extension:
            with gzip.open(genome_file_path, 'r') as genome_file:
                for chr, seq in itertools.zip_longest(*[genome_file] * 2):
                    chr = chr.decode("ascii").rstrip()[1:]
                    genome[chr] = seq.decode("ascii").rstrip()
        else:
            with open(genome_file_path, 'r') as reference_genome_file:
                with open(genome_file_path, 'r') as genome_file:
                    for chr, seq in itertools.zip_longest(*[genome_file] * 2):
                        chr = chr.rstrip()[1:]
                        genome[chr] = seq.rstrip()
        return genome

    @staticmethod
    def create_chr_ploidy_data(chr_ploidy_file_path):
        """
        Parses the chr_ploidy_data from its tab delimited resource file into a dictionary of dictionaries like so:
        {
          '1': {'male': 2, 'female': 2},
          'X': {'male': 1, 'female': 2}.
          ...
        }
        :param chr_ploidy_file_path: full path to the chr_ploidy data file
        :return: chr_ploidy_data expressed as a dictionary of dictionary as shown above.
        """
        df = pd.read_csv(chr_ploidy_file_path, sep='\t')
        return df.set_index('chr').to_dict(orient='index')

    @staticmethod
    def compare_genome_sequence_lengths(reference_file_path, genome_1_file_path, genome_2_file_path, chromosomes):
        comparison = {chromosome:[] for chromosome in chromosomes}
        genome = ExpressionUtils.create_genome(reference_file_path)
        [comparison[chromosome].append(len(sequence)) for chromosome, sequence
        in genome.items() if chromosome in chromosomes]
        genome = ExpressionUtils.create_genome(genome_1_file_path)
        for chromosome in chromosomes:
            seqeunce_length = len(genome.get(chromosome, ''))
            comparison[chromosome].append(seqeunce_length)
        genome = ExpressionUtils.create_genome(genome_2_file_path)
        for chromosome in chromosomes:
            seqeunce_length = len(genome.get(chromosome, ''))
            comparison[chromosome].append(seqeunce_length)
        return comparison


if __name__ == "__main__":
    #ExpressionUtils.create_genome("../../resources/index_files/hg38/Homo_sapiens.GRCh38.reference_genome.fa")
    chr_ploidy_data_ = ExpressionUtils.create_chr_ploidy_data('../../resources/index_files/GRCh38/Homo_sapiens.GRCh38.chr_ploidy.txt')
    reference_genome_file_path = '../../resources/index_files/GRCh38/Homo_sapiens.GRCh38.reference_genome.fa.gz'
    data_directory_path = '../../data/pipeline_results_run89/expression_pipeline/data'
    sample_data_folder = os.path.join(data_directory_path, f'sample1')
    results = ExpressionUtils.compare_genome_sequence_lengths(reference_genome_file_path,
                                                              os.path.join(sample_data_folder, 'custom_genome_1.fa'),
                                                              os.path.join(sample_data_folder, 'custom_genome_2.fa'),
                                                              chr_ploidy_data_.keys())
    df = pd.DataFrame.from_dict(results, orient='index', columns=['Reference Genome', 'Genome 1', 'Genome 2'])
    print(df)
