"""
Given a molecule file (as output by CAMPAREE), generate a coverage file
"""
import argparse
import collections
import re

parser = argparse.ArgumentParser(description="Given a molecule file (as output by CAMPAREE), generate a coverage file in the bedGraph format")
parser.add_argument("molecule_file", help="path to the input molecule file")
parser.add_argument("output_file", help="prefix for the output coverage files, two are output as output_file.forward.cov and output_file.reverse.cov")

args = parser.parse_args()


coverages = collections.defaultdict(collections.Counter)
with open(args.molecule_file) as molecule_file:
    for line in molecule_file:
        if line.startswith("#"):
            continue #Skip comment/header line

        
        transcript_id, chrom, start, cigar, ref_start, ref_cigar, strand, sequence = line.split('\t')
        # Really just need the ref_start, ref_cigar, chromosome and strand for coverage
        ref_start = int(ref_start)
        cigar_parts = re.findall("[\d]+[MINDSH]", ref_cigar)

        reference_idx = ref_start
        for part in cigar_parts:
            match_length = int(part[:-1])
            match_type = part[-1]
            if match_type == "M":
                # Count the match and move along the reference
                for i in range(match_length):
                    coverages[chrom,strand][reference_idx] += 1
                    reference_idx += 1
            elif match_type in "ND":
                # No match but we still move positions on the reference
                reference_idx += match_length
            else:
                # do nothing
                pass

with open(args.output_file + ".forward.cov", "w") as fwd_file, open(args.output_file + ".reverse.cov", "w") as rev_file:
    # Output headers
    fwd_header = f'track type=bedGraph name="Coverage  {args.molecule_file} forward strand" description="Coverage for {args.molecule_file} forward strand" visibility=full color=0,0,255 priority=20\n'
    rev_header = f'track type=bedGraph name="Coverage  {args.molecule_file} reverse strand" description="Coverage for {args.molecule_file} reverse strand" visibility=full color=255,0,0 priority=20\n'
    fwd_file.write(fwd_header)
    rev_file.write(rev_header)

    # Output the coverage of each strand of each chromosome
    # sorted alphabetically
    for chrom,strand in sorted(coverages.keys()):
        coverage = coverages[chrom,strand]
        if strand == "+":
            cov_file = fwd_file
        else:
            cov_file = rev_file

        block_start = None
        block_height = None
        for i in sorted(coverage.keys()):
            height = coverage[i]
            if block_height is None:
                # Start our first block
                block_start = i
                block_height = height
            elif block_height == height:
                # Still part of the existing block
                continue
            else:
                block_end = i - 1
                # Write out the block to the bed file
                cov_file.write('\t'.join([chrom,
                                        str(block_start),
                                        str(block_end),
                                        str(block_height)]) + "\n")

                # Start a new block
                block_start = i
                block_height = height

        if block_height is not None:
            # Output the last block
                block_end = i
                # Write out the block to the bed file
                cov_file.write('\t'.join([chrom,
                                        str(block_start),
                                        str(block_end),
                                        str(block_height)]) + "\n")
