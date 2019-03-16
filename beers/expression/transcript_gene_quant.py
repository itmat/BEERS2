import argparse
import re
import sys
import os
import collections


OUTPUT_TRANSCRIPT_FILE_NAME = "transcript_quantifications.txt"
OUTPUT_GENE_FILE_NAME = "gene_quantifications.txt"
OUTPUT_PSI_VALUE_FILE_NAME = "isoform_psi_value_quantifications.txt"

class TranscriptGeneQuantificationStep:
    """
    This class contains scripts to output quantification of transcripts and genes.
    """

    def __init__(self, geneinfo_filename, sample_directory):
        """
        The object is constructed with 2 input file sources (gene info, alignment to transcriptome).
        There are 4 output files, one each for quantification  of the genomic features: transcript, gene,
        intronic region, intergenic region.

        Based on the genomic feature being quantified, the output file is prefixed with the name of the
        aligned file and is appended with '_transcript_dist' - transcript, '_gene_dist'- gene,
        '_intron_dist' - intronic region, '_intergenic_dist' - intergenic region.
        :param geneinfo_filename: input information about the genes - fields are (chromosome, strand, start,
        end, exon count, exon starts, exon ends, gene name)
        :param align_filename: input information about the alignments to the transcriptome created using
        gene_files_preparation.py
        """

        self.geneinfo_filename = geneinfo_filename
        self.sample_directory = sample_directory


        # Create transcript distribution file and ensure that it doesn't currently exist
        self.transcript_dist_filename = os.path.join(sample_directory, OUTPUT_TRANSCRIPT_FILE_NAME)
        self.gene_dist_filename = os.path.join(sample_directory, OUTPUT_GENE_FILE_NAME)
        self.psi_value_dist_filename = os.path.join(sample_directory, OUTPUT_PSI_VALUE_FILE_NAME)


        # Dictionaries to keep track of length of transcript, number of uniquely mapped reads to transcript,
        # and final count of reads mapped to transcript
        # This procedure does not map all gene info keys used.  Consequently we need to insure that
        # assignments using new keys are initialized to 0

        self.transcript_quant_filename = os.path.join(sample_directory, f"transcriptome_1_kallisto_quant", "abundance.tsv")
        self.transcript_final_count = collections.defaultdict(float)
        self.transcript_gene_map = collections.defaultdict(str)
        self.psi_value_map = collections.defaultdict(list)


    def create_transcript_gene_map(self):
        # Create dictionary to map transcript id to gene id using geneinfo file
        with open(self.geneinfo_filename, 'r') as geneinfo_file:
            next(geneinfo_file)
            for line in geneinfo_file:
                fields = line.strip('\n').split('\t')
                self.transcript_gene_map[fields[7]] = fields[8]




    def make_transcript_gene_dist_file(self):

        with open(self.transcript_quant_filename, "r") as transcript_quant_file:
            for line in transcript_quant_file:
                if line.startswith("target_id"):
                    continue
                line = line.rstrip('\n').split('\t')
                transcript_id = line[0].split(':')[0]
                transcript_tpm = float(line[4])
                self.transcript_final_count[transcript_id] = transcript_tpm

        # Write the transcript quantification information to transcript quant filename
        with open(self.transcript_dist_filename, 'w') as transcript_dist_file:
            transcript_dist_file.write('#transcript_id' + '\t' + 'cnt' + '\n')
            for key, value in self.transcript_final_count.items():
                transcript_dist_file.write(str(key) + '\t' + str(round(value, 3)) + '\n')

        # Transcript : parent gene dictionary
        self.create_transcript_gene_map()

        # Add the transcript counts to parent gene to get gene count
        gene_count = collections.defaultdict(float)
        for key, value in self.transcript_final_count.items():
            gene_count[self.transcript_gene_map[key]] += value

        gene_count = collections.OrderedDict(sorted(gene_count.items()))
        
        # Write gene quantification information to gene quant filename
        with open(self.gene_dist_filename, 'w') as gene_dist_file:
            gene_dist_file.write('#gene_id' + '\t' + 'cnt' + '\n')
            for key, value in gene_count.items():
                gene_dist_file.write(str(key) + '\t' + str(round(value,3)) + '\n')

        # Create dictionary with key gene_id and values isoforms and their psi values
        for transcript_id in self.transcript_final_count.keys():
            gene_id = self.transcript_gene_map[transcript_id]
            if gene_count[gene_id]:
                self.psi_value_map[gene_id].append(transcript_id + ':' \
                    + str(self.transcript_final_count[transcript_id]/gene_count[gene_id]))
            else:
                self.psi_value_map[gene_id].append(transcript_id + ':' + str(0))

        # Write psi value information for each gene
        with open(self.psi_value_dist_filename, 'w') as psi_value_dist_file:
            psi_value_dist_file.write('#gene_id' + '\t' + 'isoform_psi_value' + '\n')
            for gene_id in self.psi_value_map.keys():
                psi_value_dist_file.write(str(gene_id) + '\t' + ','.join(self.psi_value_map[gene_id]) + '\n')
                 

    @staticmethod
    def is_output_valid(job_arguments):
        # TODO
        return True


    @staticmethod
    def main():
        """
        Entry point into script. Parses the argument list to obtain all the files needed and feeds them
        to the class constructor. Calls the appropriate scripts thereafter.
        """

        parser = argparse.ArgumentParser(description='Quantifier')
        parser.add_argument('-g', '--geneinfo_filename')
        parser.add_argument('-d', '--sample_directory')
        parser.add_argument('-r', '--align_filename')
        args = parser.parse_args()

        transcript_gene_quant = TranscriptGeneQuantificationStep(args.geneinfo_filename, args.sample_directory, args.align_filename)
        transcript_gene_quant.quantify_transcript()
        transcript_gene_quant.make_transcript_gene_dist_file()


if __name__ == "__main__":
    sys.exit(TranscriptGeneQuantificationStep.main())

# Example command
# python transcript_gene_quant.py -g 'geneinfo_file.txt' -r 'Aligned.out.sam' -o 'sampleA'
