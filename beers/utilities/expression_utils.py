import re
from io import StringIO
import os
import gzip
import itertools
import pandas as pd
from beers.utilities.general_utils import BeersUtilsException

class ExpressionUtils:

    # Line format definition for annotation file
    annot_output_format = '{chrom}\t{strand}\t{txStart}\t{txEnd}\t{exonCount}\t{exonStarts}\t{exonEnds}\t{transcriptID}\t{geneID}\t{geneSymbol}\t{biotype}\n'

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

        Note, this function will add a header to the output file, marked by a
        '#' character prefix.

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

            #Print annot file header (note the '#' prefix)
            output_annot_file.write("#" + ExpressionUtils.annot_output_format.replace('{', '').replace('}', ''))

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
                            ExpressionUtils.annot_output_format.format(
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
                ExpressionUtils.annot_output_format.format(
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



class NoExonsInGTF(BeersUtilsException):
    """Raised when GTF file contains no lines with "exon" in the 3rd column (feature_type)."""
    pass



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
