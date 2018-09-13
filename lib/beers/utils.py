import pandas as pd
from molecule import Molecule
import os
import pysam

class Utils:
    base_complements = {"A":"T","T":"A","G":"C","C":"G"}

    @staticmethod
    def create_complement_strand(strand):
        """
        Simple utility to provide the complement of the strand and return it in the
        5' to 3' direction.  Note that T is used rather than U even for RNA
        :param strand: RNA/DNA strand to complement.
        :return: complement strand in 5' to 3' direction
        """
        complement_strand = ''.join(Utils.base_complements[base] for base in strand)
        return complement_strand[::-1]

    @staticmethod
    def convert_log_data_into_molecules(log_filename):
        molecules = []
        log_df = pd.read_csv(log_filename)
        log_retained_df = log_df[log_df["note"] != "removed"]
        for index, row in log_retained_df.iterrows():
            molecule = Molecule(row["id"], row["sequence"], row["start"], row["cigar"])
            molecules.append(molecule)
        return molecules

    @staticmethod
    def extract_chromosome(chromosome, genome_fasta_filename):
        genome_chr_filename = os.path.splitext(genome_fasta_filename)[0] + f"_chr{chromosome}.fa"
        with open(genome_fasta_filename, 'r') as genome_ref_file,\
                open(genome_chr_filename, "w") as genome_chr_fasta_file:
            for line in genome_ref_file:
                if line.startswith(">"):
                    genome_chromosome = line[4:].rstrip()
                    if genome_chromosome != chromosome:
                        continue
                    ref_sequence = genome_ref_file.readline()
                    genome_chr_fasta_file.write(f">chr{chromosome}\n")
                    genome_chr_fasta_file.write(ref_sequence)
                    break
            else:
                print("No match found for chromosome {chromosome}.")

    @staticmethod
    def remove_cigars_with_N_from_bam_file(original_bam_filename):
        n_counter = 0
        read_counter = 0
        new_bam_filename = os.path.splitext(original_bam_filename)[0] + f"_noNs.bam"
        original_file = pysam.AlignmentFile(original_bam_filename, 'rb')
        new_file = pysam.AlignmentFile(new_bam_filename, "wb", template=original_file)
        for read in original_file.fetch('chr19'):
            read_counter += 1
            if 'N' in read.cigarstring.upper():
                print(read.cigarstring.upper())
                n_counter += 1
                continue
            new_file.write(read)
        print(f"Removed {n_counter} reads out of a total of {read_counter} reads")
        new_file.close()
        original_file.close()




if __name__ == "__main__":
    #print(Utils.create_complement_strand("AAGTGACCTAAG"))

    #Utils.convert_log_data_into_molecules("../../data/sizing_step.log")

    #Utils.extract_chromosome('19', '../../data/preBEERS/genome_mm9_edited.fa')

    Utils.remove_cigars_with_N_from_bam_file("../../data/preBEERS/Illumina.UNT_9575.Aligned.out.chr19_only.sorted.bam")
