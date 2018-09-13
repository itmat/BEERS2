import argparse
import re
import sys
import collections


def transcript_length(geneinfo_filename):
    """
    Computes total transcript length for each transcript and relates it to its ENSEMBL ID
    :param geneinfo_filename: The full transcript structure in the genome in terms of start and end locations
     for each exon
    :return: map of ENSEMBL ID to transcript lengths
    """

    # This procedure does not map all ensembl keys used.  Consequently we need to insure that
    # assignments using new keys are initialized to 0
    transcript_length_map = collections.defaultdict(int)

    # Open file for reading and gather fields in each line into an array.
    with open(geneinfo_filename, 'r') as geneinfo_file:
        for line in geneinfo_file:
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
            transcript_length_map[fields[7]] = transcript_length

    return transcript_length_map


def transcript_umap_count(reads_filename):
    """
    Counts number of reads uniquely mapped to a transcript
    :param reads_filename: The read alignment information
    :return: map of ENSEMBL ID to unique mapper read counts of the transcript
    """

    tx_umap_count = collections.defaultdict(int)

    with open(reads_filename, 'r') as reads_file:
        for line in reads_file:

            # Skip the SAM file header
            if line.startswith('@'):
                continue

            # Skip the reverse read
            reads_file.readline()

            # Gather line's field into a list
            fields = line.rstrip('\n').split('\t')

            # Ignore unaligned reads
            if fields[2] == "*":
                continue
            # Check for flag for unique mapper
            elif line.find('NH:i:1\t') >= 0:
                tx = fields[2].split(':')[0]
                tx_umap_count[tx] += 1
    return tx_umap_count



def quantify(geneinfo_filename, reads_filename, quant_filename):
    tx_final_count = collections.defaultdict(int)

    # The NH tag in the reads file (SAM file) tells us how many locations this read is aligned to.  This
    # pattern extracts that number.
    num_hits_pattern = re.compile('(NH:i:)(\d+)')

    tx_length = transcript_length(geneinfo_filename)
    tx_umap_count = transcript_umap_count(reads_filename)

    # Open file for reading and start a line counter
    with open(reads_filename, 'r') as reads_file:
        line_count = 0

        for line in reads_file:

            # Skip the SAM file header which gives the names of all contigs (contigs are the things we
            # are aligning the reads to, usually chromosomes but in this case transcripts).
            if line.startswith('@'):
                continue

            # Skip the reverse read we are counting fragments not reads.
            reads_file.readline()

            # First usable line - bump the counter and report every 100000 line.
            line_count += 1
            if line_count % 100000 == 0:
                print(f'{line_count}', file = sys.stderr)

            # Gather line's field into an array
            fields = line.rstrip('\n').split('\t')

            # This means the read did not align anywhere so we skip this read
            if fields[2] == '*':
                continue

            tx = fields[2].split(':')[0]

            # This probably means the transcript was not in our master list of all transcript models
            #  (the geneinfo filename).  So we skip it.  Really this should not happen
            #  but just in case.
            if not tx_length.get(tx):
                continue

            # Obtain the number of hits from the NH tag
            num_hits_match = re.search(num_hits_pattern, line)
            num_hits = int(num_hits_match.group(2))

            #seqid = fields[0] This appears to be unnecessary.

            # This mean that this read is a unique mapper.  It aligns only to this transcript.  So we are
            #  incrementing a counter for this transcript based on the fact that we found a unique aligning read
            #  and since it's a uniquely aligning read there's nothing more to do.
            if num_hits == 1:
                tx_final_count[tx] += 1
                continue


            # This dictionary keeps track of which ensembl ID's we've encountered that have multipmappers mapping
            #  to them.
            ensids = dict()
            fpk = dict()
            ensids[tx] = 1

            for _ in range(num_hits - 1):

                # This and the next line read the next aignment for the current read as we have multiple hits here.
                # (i.e., same read but aligning to a different transform)
                line = reads_file.readline()

                # Don't forget we skip over the reverse read we don't need it.
                reads_file.readline()

                # Parse the line's fields into a array
                multiple_hit_fields = line.rstrip().split('\t')

                multiple_hit_tx = multiple_hit_fields[2].split(':')[0]

                # Here we're basically adding to the dictionary with keys equal to the set of isoforms to which this
                #  read aligned and we know that there are at least two since this is a multimapper.
                if multiple_hit_fields[2] != "*" and tx_length[multiple_hit_tx] > 0:
                    ensids[multiple_hit_tx] = 1

            # Initialize an accumulator.
            total_fpk = 0

            # Here we're adding up the UNIQUE signals from all isoforms this read aligns to.
            #   But we can't just add up the UNIQUE signals.  We first want to normalize for the length of each
            #  isoform because we want to be able to say, for example, that two isoforms are expressed at
            #  the same level, even if one is twice as long as the other.  Because in this case we will get
            #  twice as many reads from the longer one.  But we're getting more reads because it's longer,
            #  not because it's expressed higher.  That's why we divide by length in the sum.
            for key in ensids.keys():
                fpk[key] = tx_umap_count[key]/tx_length[key]
                total_fpk += fpk[key]

            # This means that there are NO uniquey mapping reads to any of the isoforms that the current read
            #  is mapping to.  In this case we're just going to mete out the count for this read equally to all
            #  of the isoforms it aligned to.  That's why the count is incremented by 1/numKeys.
            if total_fpk == 0:
                num_keys = len(ensids.keys())
                for key in ensids.keys():
                    tx_final_count[key] += 1/num_keys

            # In this case there are uniquely mapping reads to the isforms.  In this case we're going to mete out
            #  the count for this read proportionally to how the unique mapping reads are distributed among
            #  the isoforms
            else:
                for key in ensids.keys():
                    # This percent of the signal goes to this isoform.
                    psi = fpk[key] / total_fpk

                    # Increment by psi, the psi values for this read should add to one.
                    tx_final_count[key] += psi

    ordered_tx_final_count = collections.OrderedDict(sorted(tx_final_count.items()))
    # Write the quantfication information to quant_filename
    with open(quant_filename, 'w') as quant_file:
        for key, value in ordered_tx_final_count.items():
            quant_file.write(str(key) + '\t' + str(value) + '\n')





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Quantifier')
    parser.add_argument('-g', '--geneinfo_filename')
    parser.add_argument('-r', '--reads_filename')
    parser.add_argument('-o', '--quant_filename')
    args = parser.parse_args()
    quantify(args.geneinfo_filename, args.reads_filename, args.quant_filename)

    # Example command
    #python quantify.py -g 'geneinfo_file.txt' -r '1_Aligned.out.sam' -o '1_quant.tsv'
