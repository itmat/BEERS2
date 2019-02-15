import argparse
import re
import sys
import os
import collections


OUTPUT_TRANSCRIPT_FILE_NAME = "transcript_quantifications.txt"
OUTPUT_GENE_FILE_NAME = "gene_quantifications.txt"

class TranscriptGeneQuantificationStep:
    """
    This class contains scripts to output quantification of transcripts and genes.
    """

    def __init__(self,
                geneinfo_filename,
                sample_directory,
                align_filename):
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
        self.align_filename = align_filename
        self.sample_directory = sample_directory
        # Create output director for the feature quantified files
        try:
            os.mkdir(self.sample_directory)
        except OSError:
            pass

        align_file = os.path.basename(align_filename)

        # Create transcript distribution file and ensure that it doesn't currently exist
        self.transcript_dist_filename = os.path.join(sample_directory, OUTPUT_TRANSCRIPT_FILE_NAME)
        try:
            os.remove(self.transcript_dist_filename)
        except OSError:
            pass

        # Create gene distribution file and ensure that it doesn't currently exist
        self.gene_dist_filename = os.path.join(sample_directory, OUTPUT_GENE_FILE_NAME)
        try:
            os.remove(self.gene_dist_filename)
        except OSError:
            pass


        # Dictionaries to keep track of length of transcript, number of uniquely mapped reads to transcript,
        # and final count of reads mapped to transcript
        # This procedure does not map all gene info keys used.  Consequently we need to insure that
        # assignments using new keys are initialized to 0
        self.transcript_length_map = collections.defaultdict(int)
        self.transcript_umap_count = collections.defaultdict(int)
        self.transcript_final_count = collections.defaultdict(int)
        self.transcript_gene_map = collections.defaultdict(str)


    def create_transcript_gene_map(self):
        # Create dictionary to map transcript id to gene id using geneinfo file
        # Map '*' to '*' to account for unmapped reads in align_file
        # Create entries with suffix '_1' and '_2' for each transcript

        self.transcript_gene_map['*'] = '*'

        with open(self.geneinfo_filename, 'r') as geneinfo_file:
            next(geneinfo_file)
            for line in geneinfo_file:
                fields = line.strip('\n').split('\t')
                self.transcript_gene_map[fields[7]] = fields[8]




    def transcript_length(self):
        """
        Computes total transcript length for each transcript
        and relates it to the transcript's gene info ID
        """

        # Open file for reading and gather fields in each line into an array.
        with open(self.geneinfo_filename, 'r') as geneinfo_file:
            for line in geneinfo_file:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')

                # Starting location(s) for transcript exons - banish trailing comma
                fields[5] = re.sub(r',$', '', fields[5])

                # Ending location(s) for transcript exons - banish trailing comma
                fields[6] = re.sub(r',$', '', fields[6])

                # Arrays of start and end locations for transcript exons
                start_location = fields[5].split(",")
                end_location = fields[6].split(",")

                # Sum all exon lengths together to create overall transcript length
                transcript_length = 0
                for i in range(len(start_location)):
                    transcript_length = transcript_length + int(end_location[i]) - int(start_location[i])

                # Map transcript's ENSEMBL ID to its length
                self.transcript_length_map[fields[7]] = transcript_length



    def transcript_umap(self):
        """
        Counts number of reads uniquely mapped to a transcript
        and relates it to the transcript's gene info ID
        """

        with open(self.align_filename, 'r') as align_file:
            for line in align_file:

                # Skip the SAM file header
                if line.startswith('@'):
                    continue

                # Skip the reverse read
                align_file.readline()

                # Gather line's field into a list
                fields = line.rstrip('\n').split('\t')

                # Ignore unaligned reads
                if fields[2] == "*":
                    continue
                # Check for flag for unique mapper
                elif line.find('NH:i:1\t') >= 0:
                    transcript_id = fields[2].split(':')[0]
                    self.transcript_umap_count[transcript_id] += 1


    def quantify_transcript(self):
        # The NH tag in the reads file (SAM file) tells us how many locations this read is aligned to.  This
        # pattern extracts that number.
        num_hits_pattern = re.compile('(NH:i:)(\d+)')

        self.transcript_length()
        self.transcript_umap()

        # Open file for reading and start a line counter
        with open(self.align_filename, 'r') as align_file:
            line_count = 0

            for line in align_file:

                # Skip the SAM file header which gives the names of all contigs (contigs are the things we
                # are aligning the reads to, usually chromosomes but in this case transcripts).
                if line.startswith('@'):
                    continue

                # Skip the reverse read we are counting fragments not reads.
                align_file.readline()

                # First usable line - bump the counter and report every 100000 line.
                line_count += 1
                if line_count % 100000 == 0:
                    print("{}".format(line_count), file = sys.stderr)

                # Gather line's field into an array
                fields = line.rstrip('\n').split('\t')

                # This means the read did not align anywhere so we skip this read
                if fields[2] == '*':
                    continue

                transcript_id = fields[2].split(':')[0]

                # This probably means the transcript was not in our master list of all transcript models
                #  (the geneinfo filename).  So we skip it.  Really this should not happen
                #  but just in case.
                if not self.transcript_length_map.get(transcript_id):
                    continue

                # Obtain the number of hits from the NH tag
                num_hits_match = re.search(num_hits_pattern, line)
                num_hits = int(num_hits_match.group(2))

                #seqid = fields[0] This appears to be unnecessary.

                # This mean that this read is a unique mapper.  It aligns only to this transcript.  So we are
                #  incrementing a counter for this transcript based on the fact that we found a unique aligning read
                #  and since it's a uniquely aligning read there's nothing more to do.
                if num_hits == 1:
                    self.transcript_final_count[transcript_id] += 1
                    continue


                # This dictionary keeps track of which ensembl ID's we've encountered that have multipmappers mapping
                #  to them.
                ensids = dict()
                fpk = dict()
                ensids[transcript_id] = 1

                for _ in range(num_hits - 1):

                    # This and the next line read the next aignment for the current read as we have multiple hits here.
                    # (i.e., same read but aligning to a different transform)
                    line = align_file.readline()

                    # Don't forget we skip over the reverse read we don't need it.
                    align_file.readline()

                    # Parse the line's fields into an array
                    multiple_hit_fields = line.rstrip('\n').split('\t')

                    multiple_hit_transcript = multiple_hit_fields[2].split(':')[0]


                    # Here we're basically adding to the dictionary with keys equal to the set of isoforms to which this
                    #  read aligned and we know that there are at least two since this is a multimapper.
                    if multiple_hit_fields[2] != "*" and self.transcript_length_map[multiple_hit_transcript] > 0:
                        ensids[multiple_hit_transcript] = 1

                # Initialize an accumulator.
                total_fpk = 0

                # Here we're adding up the UNIQUE signals from all isoforms this read aligns to.
                #   But we can't just add up the UNIQUE signals.  We first want to normalize for the length of each
                #  isoform because we want to be able to say, for example, that two isoforms are expressed at
                #  the same level, even if one is twice as long as the other.  Because in this case we will get
                #  twice as many reads from the longer one.  But we're getting more reads because it's longer,
                #  not because it's expressed higher.  That's why we divide by length in the sum.
                for key in ensids.keys():
                    fpk[key] = self.transcript_umap_count[key]/self.transcript_length_map[key]
                    total_fpk += fpk[key]

                # This means that there are NO uniquey mapping reads to any of the isoforms that the current read
                #  is mapping to.  In this case we're just going to mete out the count for this read equally to all
                #  of the isoforms it aligned to.  That's why the count is incremented by 1/numKeys.
                if total_fpk == 0:
                    num_keys = len(ensids.keys())
                    for key in ensids.keys():
                        self.transcript_final_count[key] += 1/num_keys

                # In this case there are uniquely mapping reads to the isforms.  In this case we're going to mete out
                #  the count for this read proportionally to how the unique mapping reads are distributed among
                #  the isoforms
                else:
                    for key in ensids.keys():
                        # This percent of the signal goes to this isoform.
                        psi = fpk[key] / total_fpk

                        # Increment by psi, the psi values for this read should add to one.
                        self.transcript_final_count[key] += psi

        self.transcript_final_count = collections.OrderedDict(sorted(self.transcript_final_count.items()))



    def make_transcript_dist_file(self):
        # Write the transcript quantification information to transcript quant filename
        with open(self.transcript_dist_filename, 'w') as transcript_dist_file:
            transcript_dist_file.write('transcript_id' + '\t' + 'cnt' + '\n')
            for key, value in self.transcript_final_count.items():
                transcript_dist_file.write(str(key) + '\t' + str(round(value, 3)) + '\n')

    def make_gene_dist_file(self):
        # Transcript : parent gene dictionary
        transcript_gene_map = collections.defaultdict(str)
        with open(self.geneinfo_filename, 'r') as geneinfo_file:
            for line in geneinfo_file:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                transcript_gene_map[fields[7]] = fields[8]

        # Add the transcript counts to parent gene to get gene count
        gene_count = collections.defaultdict(float)
        for key, value in self.transcript_final_count.items():
            gene_count[transcript_gene_map[key]] += value

        gene_count = collections.OrderedDict(sorted(gene_count.items()))
        # Write gene quantification information to gene quant filename
        with open(self.gene_dist_filename, 'w') as gene_dist_file:
            gene_dist_file.write('gene_id' + '\t' + 'cnt' + '\n')
            for key, value in gene_count.items():
                gene_dist_file.write(str(key) + '\t' + str(round(value,3)) + '\n')


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
        transcript_gene_quant.make_transcript_dist_file()
        transcript_gene_quant.make_gene_dist_file()



if __name__ == "__main__":
    sys.exit(TranscriptGeneQuantificationStep.main())

# Example command
# python transcript_gene_quant.py -g 'geneinfo_file.txt' -r 'Aligned.out.sam' -o 'sampleA'
