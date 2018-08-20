import argparse
import re
import sys
import os


class GeneFilesPreparation:
    """
    This class contains scripts to produce a fasta file for genes given the genome fasta file, a file containing
    exon locations, and a gene info file.  Additionally, any line in the gene info file related to a chromosome not
    available in the genome fasta file, is discarded in a new version of the gene info file.
    """

    def __init__(self,
                 genome_fasta_filename,
                 exons_filename,
                 geneinfo_input_filename,
                 geneinfo_output_filename,
                 genes_fasta_filename):
        """
        The object is constructed with 3 input file sources (genome fasta, exon location listing, gene info) and
        2 output file sources (relevant gene info, genes fasta).
        :param genome_fasta_filename: input genome fasta filename containing all the chromosomes of interest.  No
        line breaks are allowed within the chromosome sequence.
        :param exons_filename: input filename of exon location strings
        :param geneinfo_input_filename: input information about the genes - fields are (chromosome, strand, start,
        end, exon count, exon starts, exon ends, gene name)
        :param geneinfo_output_filename: output - contains same info as the geneinfo_input_filename but with any
        genes referencing chromosomes not in the genome fasta file removed.
        :param genes_fasta_filename: - output - contains the genes in fasta format.
        """
        self.genome_fasta_filename = genome_fasta_filename
        self.edited_genome_fasta_filename = os.path.splitext(genome_fasta_filename)[0] + "_edited.fa"
        self.exons_filename = exons_filename
        self.geneinfo_input_filename = geneinfo_input_filename
        self.geneinfo_output_filename = geneinfo_output_filename
        self.genes_fasta_filename = genes_fasta_filename

        # Provides the regex pattern for extracting components of the exon location string
        self.exon_info_pattern = re.compile('(.*):(\d+)-(\d+)')

        # Dictionaries to record whether a chromosome is available in the genome file, is available in the exon
        # file.
        self.chromosome_in_genome_file = dict()
        self.chromosome_in_exon_file = dict()

        # Since the genes fasta file will be open for appending, we need to insure that the file doesn't
        # currently exist.
        try:
            os.remove(self.genes_fasta_filename)
        except OSError:
            pass

    def prepare_gene_files(self):

        self.scrub_genome_fasta_file()

        # Open the genome fasta file for reading only.
        with open(self.edited_genome_fasta_filename, 'r') as genome_fasta_file:

            # Iterate over each chromosome in the genome fasta file.  Note that the chromosome sequence is expected on
            # only one line at this stage.
            for line in genome_fasta_file:

                # Remove the leading '>' character to leave the chromosome
                chromosome = line.lstrip('>').rstrip('\n')

                # Note that the chromosome 'chromosome' is among those listed in the genome file
                self.chromosome_in_genome_file[chromosome] = True

                # Collect the chromosome sequence from the following line
                sequence = genome_fasta_file.readline().rstrip('\n')

                # Use the chromosome and its sequence to construct a mapping of exon location string to
                # their sequences.
                exon_sequence_map = self.create_exon_sequence_map(chromosome, sequence)
                print(f"done with exons for {chromosome}")

                # Generate the gene fasta file using the genome chromosome and the related exon sequence map.
                self.make_gene_fasta_file(chromosome, exon_sequence_map)
                print(f"done with genes for {chromosome}")

        # Finally create an undated geneinfo file with any gene unrelated to the given genome chromosomes provided,
        # discarded.
        self.update_geneinfo_file()

    def scrub_genome_fasta_file(self):
        """
        Edits the genome fasta file, creating an edited version (genome fasta filename without extension + _edited.fa).
        Edits include:
        1.  Removing suplemmental information from the description line
        2.  Removing internal newlines in the sequence
        3.  Insuring all bases in sequence are represented in upper case.
        This edited file is the one used in subsequent scripts.
        """
        in_sequence = False
        with open(self.genome_fasta_filename, 'r') as genome_fasta_file, \
                open(self.edited_genome_fasta_filename, 'w') as edited_genome_fasta_file:
            for line in genome_fasta_file:
                if line.startswith('>'):
                    identifier_only = re.sub(r'[ \t].*', '', line)
                    if in_sequence:
                        edited_genome_fasta_file.write("\n")
                        in_sequence = False
                    edited_genome_fasta_file.write(identifier_only)
                else:
                    edited_genome_fasta_file.write(line.rstrip('\n').upper())
                    in_sequence = True
            edited_genome_fasta_file.write("\n")

    def update_geneinfo_file(self):
        """
        Create an update geneinfo file in which lines related to chromosomes that are not present in the
        genome fasta file are expunged.
        """
        missing_genome_chromosomes = []
        for chromosome in self.chromosome_in_exon_file.keys():
            if chromosome not in self.chromosome_in_genome_file:
                print(f"no genome sequence for {chromosome}", file=sys.stderr)
                missing_genome_chromosomes.append(chromosome)
            else:
                print(f"sequence available for {chromosome}", file=sys.stderr)
        if missing_genome_chromosomes:
            print(f"Removing the genes on the chromosomes for which no genome sequences are available.")
        with open(self.geneinfo_input_filename, 'r') as geneinfo_in, \
                open(self.geneinfo_output_filename, 'w') as geneinfo_out:
            for line in geneinfo_in:
                if not any(chromosome in line for chromosome in missing_genome_chromosomes):
                    geneinfo_out.write(line)

    def create_exon_sequence_map(self, genome_chromosome, sequence):
        """
        For the given genome chromosome and its sequence, create a dictionary of exon sequences keyed to the exon's
        location (i.e., chr:start-end).
        :param genome_chromosome:
        :param sequence:
        :return:
        """
        # Start with a empty dictionary
        exon_sequence_map = dict()

        # Open for reading only, the file containing a listing of exon locations
        with open(self.exons_filename, 'r') as exons_file:

            # Iterate over the lines in the file
            for line in exons_file:

                # Remove the newline character to be left with the exon location string only
                exon_key = line.rstrip('\n')

                # Extract the chromosome, start and end from the exon location string
                exon_info_match = re.search(self.exon_info_pattern, exon_key)
                chromosome = exon_info_match.group(1)
                exon_start = int(exon_info_match.group(2))
                exon_end = int(exon_info_match.group(3))

                # Note that the exon listing contains at least one exon on chromosome 'chromosome'
                self.chromosome_in_exon_file[chromosome] = True

                # If the chromosome on which the exon is located is the same of the genome chromosome
                # provided as a parameter, get the sequence for that exon and create a dictionary
                # entry relating the exon location string to the exon sequence.
                if chromosome == genome_chromosome:
                    exon_sequence_map[exon_key] = sequence[exon_start-1:exon_end]

        # Return the dictionary of exon location string : exon sequence for the genome chromosome
        # provided in the parameter list.
        return exon_sequence_map

    def make_gene_fasta_file(self, genome_chromosome, exon_sequence_map):

        # Open the original gene info file for reading and the genes fasta file for appending.
        with open(self.geneinfo_input_filename, 'r') as geneinfo_file, \
             open(self.genes_fasta_filename, 'a') as genes_fasta_file:

                # Iterate over the gene info file
                for line in geneinfo_file:

                    # Collect all the field values for the line read following newline removal.
                    (chromosome, strand, start, end, exon_count, exon_starts, exon_ends, name) =\
                        line.rstrip('\n').split('\t')

                    # Remove trailing commas in the exon starts and exon ends fields and split
                    # the starts and stops into their corresponding lists
                    exon_starts_list = re.sub(r'\s*,\s*$', '', exon_starts).split(",")
                    exon_ends_list = re.sub(r'\s*,\s*$', '', exon_ends).split(",")

                    # If the chromosome on which this gene is located is the same of the genome chromosome
                    # provided as a parameter, render this gene as an entry in the genes fasta file
                    if chromosome == genome_chromosome:

                        # Initialize the gene's sequence
                        gene_sequence = ""

                        # For each exon belonging to the gene, construct the exon's location string and use it as
                        # a key to obtain the actual exon sequence.  Concatenate that exon sequence to the gene
                        # sequence.  Note that the 1 added to each exon start takes into account the zero based
                        # and half-open ucsc coordinates.
                        for index in range(int(exon_count)):
                            exon_key = f'{chromosome}:{(int(exon_starts_list[index]) + 1)}-{(exon_ends_list[index])}'
                            exon_sequence = exon_sequence_map[exon_key]
                            gene_sequence += exon_sequence

                        # Remove some spurious characters from the gene name
                        gene_name = re.sub(r'::::.*', '', name)
                        # TODO determine if this substitution is still needed.
                        gene_name = re.sub(r'\([^(]+$', '', gene_name)

                        # Write the 1st line of the fasta entry - gene location string
                        genes_fasta_file.write(f'>{gene_name}:{chromosome}:{start}-{end}_{strand}\n')

                        # Write the 2nd line of the fasta entry - gene sequence.
                        genes_fasta_file.write(gene_sequence + '\n')

    @staticmethod
    def main():
        """
        Entry point into script.  Parses the argument list to obtain all the files needed and feeds them
        to the class constructor.  Finally calls the main script.
        """
        parser = argparse.ArgumentParser(description='Create Gene Info File')
        parser.add_argument('-f', '--genome_fasta_filename')
        parser.add_argument('-e', '--exons_filename')
        parser.add_argument('-i', '--geneinfo_input_filename')
        parser.add_argument('-o', '--geneinfo_output_filename')
        parser.add_argument('-g', '--genes_fasta_filename')

        args = parser.parse_args()
        file_prep = GeneFilesPreparation(args.genome_fasta_filename,
                                         args.exons_filename,
                                         args.geneinfo_input_filename,
                                         args.geneinfo_output_filename,
                                         args.genes_fasta_filename)
        file_prep.prepare_gene_files()


if __name__ == "__main__":
    sys.exit(GeneFilesPreparation.main())

'''
Usage Example:
python gene_files_preparation.py -f ../../data/genome_mm10.fa -e ../../data/exons_mm10.txt -i \
  ../../data/geneinfo_mm10.txt  -o ../../data/new_geneinfo_mm10.txt -g ../../data/genes_mm10.fa
'''
