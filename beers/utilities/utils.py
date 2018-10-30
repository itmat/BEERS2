import os
import re
import sys
import pandas as pd
import time as time
from ..library_prep.molecule import Molecule

class Utils:
    """
    Convenience methods that may find application in BEERS.
    """

    base_complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N":"N"}
    # TODO: should we allow 'N's to be complemented?

    #Line format definition for annotation file
    annot_output_format = '{chrom}\t{strand}\t{txStart}\t{txEnd}\t{exonCount}\t{exonStarts}\t{exonEnds}\t{transcriptID}\t{geneID}\t{geneSymbol}\t{biotype}\n'

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
        """
        Use BEERS pipeline log file to create molecule objects based upon the log file output as a means of running
        one step within the pipeline without having to run the entire pipeline.
        :param log_filename: path to logfile that contains molecule data that could be feed into a downstream step.
        :return: list of molecules to be used as input to a BEERS step.
        """
        molecules = []
        log_df = pd.read_csv(log_filename)

        # Do not want to pick up removed molecules
        log_retained_df = log_df[log_df["note"] != "removed"]

        # Use the fields in the log file to construct molecule objects
        for index, row in log_retained_df.iterrows():
            molecule = Molecule(row["id"], row["sequence"], row["start"], row["cigar"])
            molecules.append(molecule)
        return molecules

    @staticmethod
    def extract_chromosome(chromosome, genome_fasta_filename):
        """
        Extracts the given chromosome from the genome fasta file and outputs it into a separate fasta file.  This
        assumes that the chromosome line starts with >chr.  Further, this assume that the genome fasta file has the
        chromosome sequence all on one line.
        :param chromosome: chromosome to locate
        :param genome_fasta_filename: fasta file containing multiple chromosomes
        """

        # Create a name for the output file from the genome fasta filename by suffixing it with '_chr'.
        genome_chr_filename = os.path.splitext(genome_fasta_filename)[0] + f"_chr{chromosome}.fa"

        # Open the genome fasta file for reading and the chromosome fasta file for writing.
        with open(genome_fasta_filename, 'r') as genome_ref_file,\
                open(genome_chr_filename, "w") as genome_chr_fasta_file:

            # Iterate over the genome fasta file to find the chromosome.  Note this it is assumed that chromosome
            # lines start with '>chr'.
            for line in genome_ref_file:
                if line.startswith(">"):
                    genome_chromosome = line[4:].rstrip()
                    if genome_chromosome != chromosome:
                        continue

                    # Once the chromosome is identified, write it and its sequence into the chromosome fasta file.
                    ref_sequence = genome_ref_file.readline()
                    genome_chr_fasta_file.write(f">chr{chromosome}\n")
                    genome_chr_fasta_file.write(ref_sequence)
                    break
            else:
                print("No match found for chromosome {chromosome}.")

    @staticmethod
    def scrub_genome_fasta_file(genome_fasta_filename):
        """
        For each chromosome in the genome fasta filename provided, assemble the chromosome sequence on
        a single line and return the result as a file in the same file path as the given filename with
        'edited' suffixed to the given filename prior to the file extension.
        :param genome_fasta_filename: genome fasta file
        """

        # Create the path name for the output file.
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
    def compare_list_orders(values1, values2, missings_allowed=True):
        """
        Compare two lists to determine whether the elements that are common to both are
        ordered identically in both list.  If missingsAllowed is false, all elements must
        be common to both
        :param values1: first list to be compared
        :param values2: second list to be compared
        :param missings_allowed: whether missing elements are allowed in both lists.  Defaults to True
        :return: True if the lists are ordered identically and False otherwise.
        """

        # Empty lists or None types throw Exceptions
        if not values1:
            raise Exception(f"The first list provided is empty.")
        if not values2:
            raise Exception(f"The second list provided is empty.")

        # If missing elements are not allowed we can zip the two list together and compare every
        # element at every position.  If the lists are identically ordered the list comprehension will
        # return no discrepancies.
        if not missings_allowed:
            discrepancies = [f'value 1: {value1}, value 2: {value2}'
                             for value1, value2 in zip(values1, values2) if value1 != value2]

        # If missing element are allowed, we remove them from each list without disturbing the
        # ordering of those remaining.  If the lists containing common elements are identically
        # ordered the list comprehension will return no discrepancies.
        else:
            values1_common = [value1 for value1 in values1 if value1 in values2]
            values2_common = [value2 for value2 in values2 if value2 in values1]
            discrepancies = [(value1, value2) for value1, value2 in zip(values1_common, values2_common)
                             if value1 != value2]

        # If we have descrepancies, report them via stderr and return False.  Otherwise return True.
        if discrepancies:
            sys.stderr.write(f"{discrepancies}\n")
            return False
        else:
            return True

    @staticmethod
    def compare_file_field_orders(filename1, filename2, missings_allowed=True, delimiters=('\t', '\t'), fields=(0, 0)):
        """
        Pulls values out of a given field (defaults 0s) using the given field delimiter (default tabs) and creates
        lists of the unique values in that field for each of the two files given by their file paths and
        determines whether those lists exhibit the same order (missing elements are ignored if the
        missingAllowed parameter is set to True - the default).  The assumption is that the common values are
        grouped together on neighboring lines.  If that is found not the case, an exception is raised.
        :param filename1: file path of the first file for comparison
        :param filename2: file path of the second file for comparison
        :param missings_allowed: whether missing elements are allowed in the lists generated from both file.
        Default is True
        :param delimiters: Delimiters used by the files to separate fields.  First delimiter for first file and
        second delimiter for second file.  Defaults to tabs for both.
        :param fields: Positions of the field to use.  First position for first file and second position for
        second file.
        :return: True if the lists generated from the given files are ordered identically and False otherwise.
        """

        # Initialize the first list and iterate over the lines in the first file.
        values1 = []
        with open(filename1, 'r') as file1:
            for line in file1:

                # Grab the value in the specified field
                value = line.rstrip('\n').split(delimiters[0])[fields[0]]

                # If the value is not already in the associated list, add it.
                if value not in values1:
                    values1.append(value)

                # If the value is already in the associated list, insure that it was the last
                # value added to the list.  If not, common values are not grouped together, so
                # throw an exception.
                else:
                    if values1[len(values1)-1] != value:
                        raise Exception(f"{value} in file {filename1} is not grouped together.")

        # Initialize the second list and iterate over the lines in the second file.
        values2 = []
        with open(filename2, 'r') as file2:
            for line in file2:

                # Grab the value in the specified field
                value = line.split(delimiters[1])[fields[1]]

                # If the value is not already in the associated list, add it.
                if value not in values2:
                    values2.append(value)

                # If the value is already in the associated list, insure that it was the last
                # value added to the list.  If not, common values are not grouped together, so
                # throw an exception.
                else:
                    if values2[len(values2)-1] != value:
                        raise Exception(f"{value} in file {filename2} is not grouped together.")

        # Check whether the generated lists have identical ordering.
        return Utils.compare_list_orders(values1, values2, missings_allowed)

    @staticmethod
    def retrieve_chromosome_names_from_fasta_file(input_filename, output_filename):
        """
        Pull the chromosome names out of the given fasta file and transfer them to the given output file.
        :param input_filename: input fasta file
        :param output_filename: resulting output file containing a chromosome name on each line in the same
        order as was found in the input file.
        """
        # Pull out each line containing a chromosome name and append it to the chromosomes list.
        chromosomes = []
        with open(input_filename, 'r') as input_file:
            for line in input_file:
                if line[0] == '>':
                    chromosomes.append(line[1:].rstrip('\n'))

        # Write the list of chromosomes retrieved into the output file, one per line.
        with open(output_filename, 'w') as output_file:
            for chromosome in chromosomes:
                output_file.write(f"{chromosome}\n")

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

            #Print annot file header
            output_annot_file.write(Utils.annot_output_format.replace('{', '').replace('}', ''))

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
                            Utils.annot_output_format.format(
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
                Utils.annot_output_format.format(
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

class BeersUtilsException(Exception):
    """Base class for other Utils exceptions."""
    pass

class NoExonsInGTF(BeersUtilsException):
    """Raised when GTF file contains no lines with "exon" in the 3rd column (feature_type)."""
    pass


if __name__ == "__main__":
    # print(Utils.create_complement_strand("AAGTGACCTAAG"))

    # Utils.convert_log_data_into_molecules("../../data/sizing_step.log")

    # Utils.extract_chromosome('19', '../../data/preBEERS/genome_mm9_edited.fa')

    # Utils.remove_cigars_with_n_from_bam_file("../../data/preBEERS/Illumina.UNT_9575.Aligned.out.chr19_only.sorted.bam")

    # Utils.scrub_genome_fasta_file("../../data/expression/GRCh38/Homo_sapiens.GRCh38.reference_genome.fa")

    # list1 = [1, 5, 10, 3, 4]
    # list2 = [2, 1, 4, 11, 7, 6]
    # print(Utils.compare_list_orders(list1, list2))

    # Utils.retrieve_chromosome_names_from_fasta_file('../../data/preBEERS/hg19/reference_genome_edited.fa',
    #                                                 '../../data/preBEERS/hg19/reference_chromosomes.txt')
    #print(Utils.compare_file_field_orders('../../data/preBEERS/hg19/reference_chromosomes.txt',
    #                                      '../../data/preBEERS/hg19/Test_dataset.All_unique_mappers'
    #                                      '.fw_only_variants.txt',
    #                                      delimiters=(' ', ':')))

    #print(Utils.generate_seed())

    pass
