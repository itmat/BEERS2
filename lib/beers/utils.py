import pandas as pd
from molecule import Molecule
import os
import pysam
import re

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

    @staticmethod
    def scrub_genome_fasta_file(genome_fasta_filename):

        edited_genome_fasta_filename = os.path.splitext(genome_fasta_filename)[0] + "_edited.fa"

        in_sequence = False

        with open(genome_fasta_filename, 'r') as genome_fasta_file, \
                open(edited_genome_fasta_filename, 'w') as edited_genome_fasta_file:

            # Iterate over the lines in the original genome fasta file
            for line in genome_fasta_file:

                # Identify whether the current line is a description line or a sequence line
                if line.startswith('>'):

                    # For a description line, remove any supplemental information following the identifier and
                    # if the in_sequence flag is raised, lower it and add a line break to the new genome fasta
                    # file before adding the modified description line.
                    identifier_only = re.sub(r'[ \t].*', '', line)
                    if in_sequence:
                        edited_genome_fasta_file.write("\n")
                        in_sequence = False
                    edited_genome_fasta_file.write(identifier_only)
                # Otherwise, add the sequence to the new genome fasta file after removing the line break and
                # insuring all bases are in upper case.  Also raise the in sequence flag.
                else:
                    edited_genome_fasta_file.write(line.rstrip('\n').upper())
                    in_sequence = True

            # Finally add a line break to the end of the new genome fasta file.
            edited_genome_fasta_file.write("\n")

    @staticmethod
    def convert_gtf_to_annot_file_format(gtf_filename):
        """Convert a GTF file to a tab-delimited annotation file with one line
        per transcript. Each line in the annotation file will have the following
        columns:
             1 - chrom
             2 - strand
             3 - txStart
             4 - txEnd
             5 - exonCount
             6 - exonStarts
             7 - exonEnds
             8 - transcript_id
             9 - gene_id
            10 - genesymbol
            11 - biotype
        This method derives transcript info from the "exon" lines in the GTF
        file and assumes the exons are listed in the order they appear in the
        transcript, as opposed to their genomic coordinates. The annotation file
        will list exons in order by plus-strand coordinates, so this method
        reverses the order of exons for all minus-strand transcripts.

        See website for standard 9-column GTF specification:
        https://useast.ensembl.org/info/website/upload/gff.html

        Parameters
        ----------
        gtf_filename : string
            Path to GTF file to be converted annotation file format.

        Returns
        -------
        string
            Name of the annotation file produced from the GTF file.

        """
        #Note, as-is this statement won't expand "~" in the path to point to
        #home directory.
        output_annot_filename = os.path.splitext(gtf_filename)[0] + ".annotation.txt"

        with open(gtf_filename, 'r') as gtf_file, \
                open(output_annot_filename, 'w') as output_annot_file:

            #Define line format for annotation file
            annot_output_format = '{chrom}\t{strand}\t{txStart}\t{txEnd}\t{exonCount}\t{exonStarts}\t{exonEnds}\t{transcriptID}\t{geneID}\t{geneSymbol}\t{biotype}\n'

            #Print annot file header
            output_annot_file.write(annot_output_format.replace('{', '').replace('}', ''))

            #Regex patterns used to extract individual attributes from the 9th
            #column in the GTF file (the "attributes" column)
            txid_pattern = re.compile(r'transcript_id "([^"]+)";')
            geneid_pattern = re.compile(r'gene_id "([^"]+)";')
            genesymbol_pattern = re.compile(r'gene_name "([^"]+)";')
            biotype_pattern = re.compile(r'gene_biotype "([^"]+)";')

            line_data = [] #List of line fields from current line of gtf
            curr_gtf_tx = "" #transcript ID from current line of gtf
            chrom = "" #Chromosome for current transcript
            strand = "" #Strand for current transcript
            txid = "" #ID of current transcript
            geneid = "" #gene ID of current transcript
            genesymbol = "None" #gene symbol of current transcript
            biotype = "None" #Biotype of current transcript
            ex_starts = [] #List of exon start coordinates for current transcript
            ex_stops = [] #List of exon stop coordinates for current transcript
            ex_count = 1 #Number of exons in current transcript

            #Step through GTF until first exon entry (need to prime variables
            #with data from first exon).
            exon_found = False
            for line in gtf_file:
                line_data = line.split("\t")
                #First conditional is to skip any header lines, which would
                #result in a "list index out of range" error when checking the
                #second conditional. Header lines don't contain any tabs, so
                #there is no 3rd list element in line_data, following the split
                #command.
                if len(line_data) > 1 and line_data[2] == "exon":
                    exon_found = True
                    break

            if not exon_found:
                raise NoExonsInGTF('ERROR: {gtf_file} contains no lines with exon feature_type.\n'.format(gtf_file=gtf_filename))

            #Prime variables with data from first exon
            chrom = line_data[0]
            strand = line_data[6]
            ex_starts.append(line_data[3])
            ex_stops.append(line_data[4])
            txid = txid_pattern.search(line_data[8]).group(1)
            geneid = geneid_pattern.search(line_data[8]).group(1)
            if genesymbol_pattern.search(line_data[8]):
                genesymbol = genesymbol_pattern.search(line_data[8]).group(1)
            if biotype_pattern.search(line_data[8]):
                biotype = biotype_pattern.search(line_data[8]).group(1)

            #process the remainder of the GTF file
            for line in gtf_file:
                line = line.rstrip('\n')
                line_data = line.split("\t")

                if line_data[2] == "exon":
                    curr_gtf_tx = txid_pattern.search(line_data[8]).group(1)

                    #Check transcript in current line is a new transcript
                    if curr_gtf_tx != txid:

                        #Check strand of transcript. If minus, reverse exon order.
                        if strand == "-":
                            ex_starts.reverse()
                            ex_stops.reverse()

                        #Format data from previous transcript and write to annotation file
                        output_annot_file.write(
                            annot_output_format.format(
                                chrom=chrom,
                                strand=strand,
                                txStart=ex_starts[0],
                                txEnd=ex_stops[-1],
                                exonCount=ex_count,
                                exonStarts=','.join(ex_starts),
                                exonEnds=','.join(ex_stops),
                                transcriptID=txid,
                                geneID=geneid,
                                geneSymbol=genesymbol,
                                biotype=biotype
                            )
                        )

                        #Load data from new transcript into appropriate variables
                        txid = curr_gtf_tx
                        chrom = line_data[0]
                        strand = line_data[6]
                        ex_count = 1
                        ex_starts = [line_data[3]]
                        ex_stops = [line_data[4]]
                        geneid = geneid_pattern.search(line_data[8]).group(1)
                        if genesymbol_pattern.search(line_data[8]):
                            genesymbol = genesymbol_pattern.search(line_data[8]).group(1)
                        else:
                            genesymbol = "None"
                        if biotype_pattern.search(line_data[8]):
                            biotype = biotype_pattern.search(line_data[8]).group(1)
                        else:
                            biotype = "None"

                    #This exon is strill from the same transcript.
                    else:
                        ex_starts.append(line_data[3])
                        ex_stops.append(line_data[4])
                        ex_count += 1

            #Finish processing last transcript in GTF file

            #Check strand. If minus, reverse exon order.
            if strand == "-":
                ex_starts.reverse()
                ex_stops.reverse()

            #Format data from last transcript and write to annotation file
            output_annot_file.write(
                annot_output_format.format(
                    chrom=chrom,
                    strand=strand,
                    txStart=ex_starts[0],
                    txEnd=ex_stops[-1],
                    exonCount=ex_count,
                    exonStarts=','.join(ex_starts),
                    exonEnds=','.join(ex_stops),
                    transcriptID=txid,
                    geneID=geneid,
                    geneSymbol=genesymbol,
                    biotype=biotype
                )
            )

        return output_annot_filename

class BeersUtilsException(Exception):
    """Base class for other Utils exceptions."""
    pass

class NoExonsInGTF(BeersUtilsException):
    """Raised when GTF file contains no lines with "exon" in the 3rd column (feature_type)."""
    pass


if __name__ == "__main__":
    #print(Utils.create_complement_strand("AAGTGACCTAAG"))

    #Utils.convert_log_data_into_molecules("../../data/sizing_step.log")

    #Utils.extract_chromosome('19', '../../data/preBEERS/genome_mm9_edited.fa')

    #Utils.remove_cigars_with_N_from_bam_file("../../data/preBEERS/Illumina.UNT_9575.Aligned.out.chr19_only.sorted.bam")

    Utils.scrub_genome_fasta_file("../../data/preBEERS/hg19_chr21_22_ref.fa")
