import argparse
import re
import sys
import collections


def quantify(lengths_filename, geneinfo_filename, counts_filename, reads_filename):
    id_mapping = dict()
    lengths = collections.defaultdict(int)
    ucount = collections.defaultdict(int)
    final_count = collections.defaultdict(int)
    fpk = dict()
    num_hits_pattern = re.compile('(NH:i:)(\d+)')
    with open(lengths_filename, 'r') as lengths_file:
        for line in lengths_file:
            fields = line.rstrip('\n').split('\t')
            id_mapping[fields[0]] = fields[1]
    with open(geneinfo_filename, 'r') as geneinfo_file:
        for line in geneinfo_file:
            fields = line.rstrip('\n').split('\t')
            len_ = 0
            fields[5] = re.sub(r',$','', fields[5])
            fields[6] = re.sub(r',$','', fields[6])
            s = fields[5].split(",")
            e = fields[6].split(",")
            for i in range(len(s)):
                len_ = len_ + int(e[i]) - int(s[i])
            ensid_value = id_mapping[fields[7]]
            lengths[ensid_value] = len_
    with open(counts_filename, 'r') as counts_file:
        for line in counts_file:
            fields = line.rstrip('\n').split('\t')
            ucount[fields[0]] = int(fields[1])
    with open(reads_filename, 'r') as reads_file:
        line_count = 0
        for line in reads_file:
            if line.startswith('@'):
                continue
            reads_file.readline()
            line_count += 1
            if line_count % 100000 == 0:
                print(f'{line_count}', file = sys.stderr)
            fields = line.rstrip('\n').split('\t')
            if fields[2] == '*':
                continue
            if not lengths.get(fields[2]):
                continue
            num_hits_match = re.search(num_hits_pattern, line)
            num_hits = int(num_hits_match.group(2))
            #seqid = fields[0]
            if num_hits == 1:
                final_count[fields[2]] += 1
                continue
            ensids = {}
            ensids[fields[2]] = 1
            if num_hits > 1:
                for _ in range(num_hits - 1):

                    # Read line
                    line = reads_file.readline()

                    # Skip next line
                    reads_file.readline()

                    # Parse line fields
                    new_fields = line.rstrip().split('\t')


                    #if new_fields[2] != "*" and lengths[new_fields[2]] > 0:
                        #ensids[new_fields[2]] = 1
            total_fpk = 0
            for key in ensids.keys():
                fpk[key] = ucount[key]/lengths[key]
                total_fpk += fpk[key]
            if total_fpk == 0:
                num_keys = len(ensids.keys())
                for key in ensids.keys():
                    final_count[key] += 1/num_keys
            else:
                for key in ensids.keys():
                    psi = fpk[key] / total_fpk
                    final_count[key] += psi

    ordered_final_count = collections.OrderedDict(sorted(final_count.items()))
    for key, value in ordered_final_count.items():
        print(f'{key}\t{int(value)}')













if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Quantifier')
    parser.add_argument('-l', '--lengths_filename')
    parser.add_argument('-g', '--geneinfo_filename')
    parser.add_argument('-c', '--counts_filename')
    parser.add_argument('-r', '--reads_filename')
    args = parser.parse_args()
    quantify(args.lengths_filename, args.geneinfo_filename, args.counts_filename, args.reads_filename)

    #python quantify.py -l '/Users/crislawrence/Desktop/cris/BEERS2ENS_SORTED.txt' -g '/Users/crislawrence/Desktop/cris/simulator_config_geneinfo_mm10_ensembl_r75' -c '/Users/crislawrence/Desktop/cris/STAR_umap_counts.Liv1_4.tsv' -r '/Users/crislawrence/Desktop/cris/sim_reads_Liv1_4_Aligned.out.sam'
