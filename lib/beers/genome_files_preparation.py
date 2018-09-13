import sys
import argparse
import numpy as np
import os
import re
import pysam
import vcf
from io import StringIO
from timeit import default_timer as timer


class GenomeFilesPreparation:

    def __init__(self, genome_ref_filename, snp_calls_filename, genome_filename, log_filename, random_seed):
        np.random.seed(random_seed)
        self.genome_ref_filename = genome_ref_filename
        self.snp_calls_filename = snp_calls_filename
        self.genome_filename = genome_filename
        self.log_filename = log_filename
        self.allele_freqs_pattern = re.compile('AF=([,\d]+);')
        # TODO make force choice an input param
        #pysam.tabix_index(snp_calls_filename,preset='vcf',force=False,zerobased=True)
        try:
            os.remove(self.genome_filename)
        except OSError:
            pass

    def prepare_genome_files(self):
        with open(self.genome_ref_filename, 'r') as genome_ref_file, open(self.log_filename, 'w') as log_file:
            for line in genome_ref_file:
                if line.startswith(">"):
                    # Skipping >chr - seems fragile
                    chromosome = line[1:].rstrip()
                    if chromosome != 'chr19':
                        break
                    ref_sequence = genome_ref_file.readline().rstrip()
                    print(f"Length of reference sequence is {len(ref_sequence)}")
                    last_alt_accepted = False
                    offset = 0
                    where_ref_seq_resumes = 0
                    snp_call_counter = 0
                    vcf_reader = vcf.Reader(open(self.snp_calls_filename, 'r'))
                    sequence_str = StringIO()
                    sequence_length = 0
                    block_start = timer()
                    for snp_call in vcf_reader:
                        print(f"snp_call.CHROM is {snp_call.CHROM} and snp_call.ALT is {snp_call.ALT}")
                        if snp_call.CHROM == chromosome and snp_call.ALT and snp_call.ALT[0]:
                            snp_call_counter += 1
                            # making pos is 0 based
                            pos = snp_call.POS - 1
                            ref = snp_call.REF
                            alt = [item.sequence for item in snp_call.ALT]
                            assert ref == ref_sequence[pos:pos + len(ref)], \
                                    f"Reference genome base '{ref_sequence[pos:pos + len(ref)]}' and snp ref '{ref}' don't match" \
                                    f" for chromosome '{chromosome}' and 0 based position '{pos}'"

                            if snp_call_counter % 10000 == 0:
                                block_end = timer()
                                print(f"Up to snp_call_counter {snp_call_counter} at position {pos} at time {block_end - block_start}")
                                block_start = timer()

                            if pos < where_ref_seq_resumes and last_alt_accepted:
                                continue
                            # We need to back up to accommodate an overlapping change since the prior ALT was not used.
                            if pos < where_ref_seq_resumes and not last_alt_accepted:
                                sequence_str.truncate(pos+offset)
                                sequence_length = pos + offset
                                where_ref_seq_resumes = pos
                            choices = alt + [ref]

                            # Accounting for inaccuracy in floats
                            freqs = snp_call.aaf + [1 - min(1, sum(snp_call.aaf))]
                            selection = np.random.choice(choices, p=freqs)

                            #sequence_length += len(ref_sequence[where_ref_seq_resumes:pos]) + len(selection)
                            sequence_length += pos - where_ref_seq_resumes + len(selection)
                            sequence_str.write(ref_sequence[where_ref_seq_resumes:pos])
                            sequence_str.write(selection)

                            offset += len(selection) - len(ref)
                            log_file.write(f"Zero based position: {pos}\n")
                            log_file.write(f"\tRef seq length after ref {len(ref_sequence[:pos+len(ref)])} - selection is {selection} whereas ref was {ref}\n")
                            log_file.write(f"\tNew seq length after selection {sequence_length} and offset is {offset}\n\n")
                            assert sequence_length == offset + pos + len(ref), f"seqs out of sync at pos {pos}"
                            #assert sequence_length  == len(ref_sequence[:pos+len(ref)]) + offset, f"seqs out of sync at pos {pos}"
                            #input('?')
                            where_ref_seq_resumes = pos + len(ref)
                            last_alt_accepted = selection != ref
                        elif snp_call.CHROM != chromosome:
                            sequence = sequence_str.read() + ref_sequence[where_ref_seq_resumes:]
                            sequence_str.close()
                            break
                    with open(self.genome_filename, 'a') as genome_file:
                        genome_file.write(line)
                        genome_file.write(sequence + "\n")
                    print(f"SNP calls: {snp_call_counter}")
                    break

    def combine_snp_calls(self, vcf_reader, log_file):
        current_chrom = None
        current_pos = None
        current_ref = ''
        current_alt = []
        current_freqs = []
        for snp_call in vcf_reader:
            next_pos = snp_call.POS - 1
            #print(f"current_pos {current_pos} and next_pos {next_pos}")
            #input("continue")
            if current_pos and next_pos != current_pos:
                current_data = {'chrom': current_chrom, 'pos': current_pos, 'ref': current_ref, 'alt': current_alt, 'freqs': current_freqs}
                current_chrom = snp_call.CHROM
                current_pos = next_pos
                current_ref = snp_call.REF
                if snp_call.ALT and snp_call.ALT[0]:
                    current_alt = [item.sequence for item in snp_call.ALT]
                    current_freqs = snp_call.aaf
                else:
                    current_alt = []
                    current_freqs = []
                yield current_data
            else:
                if snp_call.ALT and snp_call.ALT[0]:
                    log_file.write(f"+++ position {current_pos} has one than 1 line - combining +++\n\n")
                    current_chrom = snp_call.CHROM
                    current_ref = snp_call.REF
                    current_pos = snp_call.POS - 1
                    current_alt.extend([item.sequence for item in snp_call.ALT])
                    current_freqs.extend(snp_call.aaf)
        return None






    @staticmethod
    def main():
        parser = argparse.ArgumentParser(description='Create Genome File')
        parser.add_argument('-r', '--genome_ref_filename')
        parser.add_argument('-c', '--snp_calls_filename')
        parser.add_argument('-g', '--genome_filename')
        parser.add_argument('-s', '--random_seed', default=None, type=int)
        parser.add_argument('-l', '--log_filename')

        args = parser.parse_args()
        file_prep = GenomeFilesPreparation(args.genome_ref_filename,
                                         args.snp_calls_filename,
                                         args.genome_filename,
                                         args.log_filename,
                                         args.random_seed)
        start = timer()
        file_prep.prepare_genome_files()
        end = timer()
        print(f"Genome Seq: {end - start}")


if __name__ == '__main__':
    sys.exit(GenomeFilesPreparation.main())

# Sample command
# python genome_files_preparation.py -r ../../data/preBEERS/genome_mm10_edited.fa -c ../../data/preBEERS/mgp.v5.merged.indels.dbSNP142.normed.vcf \
# -l ../../data/preBEERS/genome_chr1.log -g ../../data/preBEERS/genome_chr1.fa -s 100


# Sample command
# python genome_files_preparation.py -r ../../data/preBEERS/genome_mm9_edited_chr19.fa -c ../../data/preBEERS/Illumina.UNT_9575.Aligned.out.chr19_only.vcf \
# -l ../../data/preBEERS/genome_chr19.log -g ../../data/preBEERS/genome_chr19.fa -s 100
