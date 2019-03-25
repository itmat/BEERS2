import argparse
import re
import sys
import os
import collections


from pysam import AlignmentFile



OUTPUT_ALLELIC_IMBALANCE_FILE_NAME = "allelic_imbalance_quantifications.txt"

class AllelicImbalanceQuantificationStep:
    """
    This class contains scripts to output quantification of allelic imbalance
    """

    def __init__(self,
                sample_directory):
        """
        The object is constructed with
        (i) an input file source for gene info
        (ii) Root of the aligned filenames (alignment to transcriptome of each parent with suffixes '_1','_2'.)
        There is one output file with quantification information on the allelic imbalance of genes.

        The output file is called 'allelic_imbalance_quantifications.txt'
        :param geneinfo_filename: input information about the genes - fields are (chromosome, strand, start,
        end, exon count, exon starts, exon ends, gene name)
        :param align_filename_root: input information about the alignments to the two transcriptomes, one for each parent,
        created using gene_files_preparation.py
        """

        self.sample_directory = sample_directory
        self.geneinfo_filename_1 = os.path.join(self.sample_directory, 'updated_annotation_1.txt')
        self.geneinfo_filename_2 = os.path.join(self.sample_directory, 'updated_annotation_2.txt')
        self.genome_alignment_file = os.path.join(self.sample_directory, "genome_alignment.Aligned.sortedByCoord.out.bam")

        self.align_filename_1 = os.path.join(self.sample_directory, '1_Aligned.out.sam')
        self.align_filename_2 = os.path.join(self.sample_directory, '2_Aligned.out.sam')

        # Create allelic imbalance distribution file and ensure that it doesn't currently exist
        self.allele_imbalance_dist_filename = os.path.join(sample_directory, OUTPUT_ALLELIC_IMBALANCE_FILE_NAME)
        try:
            os.remove(self.allele_imbalance_dist_filename)
        except OSError:
            pass


        # Dictionaries to map transcripts to genes and keep track of final count of reads mapped to a gene
        # This procedure does not map all gene info keys used.  Consequently we need to insure that
        # assignments using new keys are initialized to 0
        self.transcript_gene_map = collections.defaultdict(str)
        self.gene_final_count = collections.defaultdict(lambda: collections.defaultdict(int))
        self.exclusive_genes = []

    
    def create_transcript_gene_map(self):
        """
        Create dictionary to map transcript id to gene id using geneinfo file
        Map '*' to '*' to account for unmapped reads in align_file
        Create entries with suffix '_1' and '_2' for each transcript
        """

        self.transcript_gene_map['*'] = '*'

        with open(self.geneinfo_filename_1, 'r') as geneinfo_file:
            next(geneinfo_file)
            for line in geneinfo_file:
                fields = line.strip('\n').split('\t')
                self.transcript_gene_map[fields[7]] = fields[8]

    def reads_to_ignore(self):
        reads_to_ignore = []
        bamfile = AlignmentFile(self.genome_alignment_file, "rb")
        num_hits_pattern = re.compile('(NH:i:)(\d+)')

        for read in bamfile.fetch(until_eof=True):
            num_hits = dict(read.tags)['NH']
            if num_hits > 1:
                reads_to_ignore.append(read.query_name)

        return reads_to_ignore


    def read_info(self, in_align_filename):
        """
        Create dictionary which maps a read id in SAM file to a dictionary with two keys 'transcript_id' and 'NM'.
        The value associated with 'transcript_id' is a list of all transcripts the read aligned to.
        The value associated with 'NM' is the corresponding edit distance information for each alignment.
        For non-mappers the transcript_id is '*' and edit distance is 100 (Make it read length).
        """

        read_info_map = collections.defaultdict(dict)


        # The NM tag in the SAM file tells us the edit distance for the alignment.
        # This pattern extracts that number.
        num_mismatches_pattern = re.compile('(NM:i:)(\d+)')

        with open(in_align_filename, 'r') as infile:
            for line in infile:
                if line.startswith('@'):
                    continue

                # read forward and reverse read
                forward = line
                reverse = next(infile)

                # Parse the fields for the forward read into an array
                fwd_fields = forward.rstrip('\n').split('\t')

                # Parse the fields for the reverse read into an array
                rev_fields = reverse.rstrip('\n').split('\t')

                fwd_transcript_id = fwd_fields[2].split(':')[0]
                rev_transcript_id = rev_fields[2].split(':')[0]

                # This means both forward and reverse reads are non-mappers
                # So store 'transcript_id' as '*' and 'NM' as 2*read_length
                if fwd_transcript_id == '*' and rev_transcript_id == '*':
                    read_info_map[fwd_fields[0]]['transcript_id'] =  '*' 
                    read_info_map[fwd_fields[0]]['NM'] =  200 
                    continue
                # Get transcript_id for mapped reads
                elif fwd_transcript_id == rev_transcript_id:
                    transcript_id = fwd_transcript_id
                else:
                    transcript_id = (fwd_transcript_id + rev_transcript_id).replace('*','')

                # This probably means the transcript was not in our master list of all transcript models
                #  (the geneinfo filename).  So we skip it.  Really this should not happen
                #  but just in case.
                if not self.transcript_gene_map.get(transcript_id):
                    continue

                # Obtain the edit distance information for the forward read
                fwd_NM_match = re.search(num_mismatches_pattern, forward)
                rev_NM_match = re.search(num_mismatches_pattern, reverse)
                if fwd_NM_match and rev_NM_match:
                    fwd_NM_count = int(fwd_NM_match.group(2))
                    rev_NM_count = int(rev_NM_match.group(2))
                    NM_count = fwd_NM_count + rev_NM_count
                elif not (fwd_NM_match and rev_NM_match):
                    NM_count = 200
                else:
                    NM_count = 100

                # Update read_info dictionary with transcript_id and corresponding edit distance
                read_info_map[fwd_fields[0]]['transcript_id'] = transcript_id 
                read_info_map[fwd_fields[0]]['NM'] = NM_count

        return read_info_map


    def quantify_allelic_imbalance(self):
        """
        This is the main step which quantifies allelic imbalance for all genes in the annotation based on
        the aligned files for parents 1 and 2.
        """

        self.create_transcript_gene_map()

        # Create read info dictionaries for each parent
        read_info_1 = self.read_info(self.align_filename_1)
        read_info_2 = self.read_info(self.align_filename_2)

        reads_to_ignore = self.reads_to_ignore()

        read_ids_1 = set(read_info_1.keys()).difference(reads_to_ignore)
        read_ids_2 = set(read_info_2.keys()).difference(reads_to_ignore)
        read_ids = read_ids_1.intersection(read_ids_2)
        read_ids_1_u = read_ids_1.difference(read_ids)
        read_ids_2_u = read_ids_2.difference(read_ids)

        for read in read_ids_1_u:
            transcript = read_info_1[read]['transcript_id']
            gene = self.transcript_gene_map[transcript]
            self.gene_final_count[gene]['1'] += 1

        for read in read_ids_2_u:
            transcript = read_info_2[read]['transcript_id']
            gene = self.transcript_gene_map[transcript]
            self.gene_final_count[gene]['2'] += 1

        for read in read_ids:
            # Transcripts to which the read mapped for each parent
            transcript_1 = read_info_1[read]['transcript_id']
            transcript_2 = read_info_2[read]['transcript_id']

            # The read did not map to any transcript in either parent
            if transcript_1 == '*' and transcript_2 == '*':
                continue
            # The read mapped to atleast one transcript in each parent
            elif transcript_1 != '*' and transcript_2 != '*':
                # Get the genes in parent 1 to which the read mapped
                gene_1 = self.transcript_gene_map[transcript_1] 
                NM_count_1 = read_info_1[read]['NM']

                # Get the genes in parent 2 to which the read mapped
                gene_2 = self.transcript_gene_map[transcript_2] 
                NM_count_2 = read_info_2[read]['NM']

                # Amongst the genes to which the read mapped,
                # there is exactly one gene in common between parent 1 and 2.
                if gene_1 == gene_2:
                    # Minimum edit distance for the mapping to the gene is the same in
                    # parent 1 and parent 2. So increment counts of both alleles of the genes by 0.5
                    if NM_count_1 == NM_count_2:
                        self.gene_final_count[gene_1]['1'] += 0.5
                        self.gene_final_count[gene_1]['2'] += 0.5
                    # Minimum edit distance for the mapping to the gene is less in parent 1.
                    # So increment count of allele of gene corresponding to parent 1.
                    elif NM_count_1 < NM_count_2:
                        self.gene_final_count[gene_1]['1'] += 1
                    # Minimum edit distance for the mapping to the gene is less in parent 2.
                    # So increment count of allele of gene corresponding to parent 2.
                    else:
                        self.gene_final_count[gene_1]['2'] += 1
            # The read is a non-mapper for the parent 1 transcriptome
            elif transcript_1 == '*':
                # Get the genes in parent 2 to which the read mapped
                gene_2 = self.transcript_gene_map[transcript_2]
                self.gene_final_count[gene_2]['2'] += 1

            # The read is a non-mapper for the parent 2 transcriptome
            elif transcript_2 == '*':
                gene_1 = self.transcript_gene_map[transcript_1]
                self.gene_final_count[gene_1]['1'] += 1

        #self.gene_final_count = collections.OrderedDict(sorted(self.gene_final_count.items()))


    def make_allele_imbalance_dist_file(self):
        genelist_1 = []
        with open(self.geneinfo_filename_1, 'r') as geneinfo_file_1:
            for line in geneinfo_file_1:
                if line.startswith('#'):
                    continue
                fields = line.strip('\n').split('\t')
                genelist_1.append(fields[8])

        genelist_2 = []
        with open(self.geneinfo_filename_2, 'r') as geneinfo_file_2:
            for line in geneinfo_file_2:
                if line.startswith('#'):
                    continue
                fields = line.strip('\n').split('\t')
                genelist_2.append(fields[8])

        exclusive_genes = list(set(genelist_1).difference(set(genelist_2)))

        # Write the allelic imbalance quantification information to allele imbalance dist filename
        with open(self.allele_imbalance_dist_filename, 'w') as allele_imbalance_dist_file:
            allele_imbalance_dist_file.write('#gene_id' + '\t' + '_1' + '\t' + '_2' + '\n')

            #for key, value in list(self.gene_final_count.items()):
            for gene_id in sorted(set(self.transcript_gene_map.values())):
                if gene_id in exclusive_genes:
                    allele_imbalance_dist_file.write(str(gene_id) + '\t' + str(1.0) + '\t' + str(0.0) + '\n')
                    continue

                if gene_id == "*":
                    continue

                read_count_1 = self.gene_final_count[gene_id]['1']
                read_count_2 = self.gene_final_count[gene_id]['2']
                gene_read_count = read_count_1 + read_count_2

                if gene_read_count == 0:
                    allele_imbalance_dist_file.write(str(gene_id) + '\t' + str(0.5) + '\t' + str(0.5) + '\n')
                else:
                    allele_imbalance_dist_file.write(str(gene_id) + '\t' + str(read_count_1/gene_read_count) + '\t' +\
                        str(read_count_2/gene_read_count) + '\n')



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
        parser.add_argument('-d', '--sample_directory')
        args = parser.parse_args()

        allelic_imbalance_quant = AllelicImbalanceQuantificationStep(args.sample_directory)
        allelic_imbalance_quant.quantify_allelic_imbalance()
        allelic_imbalance_quant.make_allele_imbalance_dist_file()



if __name__ == "__main__":
    sys.exit(AllelicImbalanceQuantificationStep.main())

# Example command
# python allelic_imbalance_quant.py -g 'geneinfo_file.txt' -d 'sampleA' -r 'Aligned.out'


