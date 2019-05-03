"""
Given a molecule file (as output by CAMPAREE), generate a coverage file
"""
import argparse
import collections
import re

import numpy

parser = argparse.ArgumentParser(description="Given a molecule file (as output by CAMPAREE), generate a coverage file in the bedGraph format")
parser.add_argument("molecule_file", help="path to the input molecule file")
parser.add_argument("output_file", help="prefix for the output coverage files, two are output as output_file.forward.cov and output_file.reverse.cov")
parser.add_argument("-m", "--max_chromosome_size", help="size of largest chromosome in megabases (250 for humans, 196 for mouse), this can be an overestimate, determines memory use", default=260, type=int)

args = parser.parse_args()

cigar_re = re.compile("[\d]+[MINDSH]")

def make_chrom():
    print(f"Creating the {len(coverages)+1} chromosome")
    return numpy.zeros(args.max_chromosome_size * 1_000_000, dtype="int32")
coverages = collections.defaultdict(make_chrom)
with open(args.molecule_file) as molecule_file:
    for line in molecule_file:
        if line.startswith("#"):
            continue #Skip comment/header line


        transcript_id, chrom, start, cigar, ref_start, ref_cigar, strand, sequence = line.split('\t')
        # Really just need the ref_start, ref_cigar, chromosome and strand for coverage
        ref_start = int(ref_start)
        cigar_parts = cigar_re.findall(ref_cigar)

        reference_idx = ref_start
        for part in cigar_parts:
            match_length = int(part[:-1])
            match_type = part[-1]
            if match_type == "M":
                # Count the match and move along the reference
                coverages[chrom,strand][reference_idx:reference_idx+match_length] += 1
                reference_idx += match_length
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
        print(f"Starting chrom {chrom} strand {strand}")
        coverage = coverages[chrom,strand]
        # We allocated too much space for every chromosome and most didn't use it
        # unused is left as zeros, so we can just drop those on the end
        coverage = numpy.trim_zeros(coverage, 'b')
        print(f"\t{len(coverage)} bases to output")

        if strand == "+":
            cov_file = fwd_file
        else:
            cov_file = rev_file

        block_start = None
        block_height = None
        for i, height in enumerate(coverage):
            if block_height is None:
                # Start our first block
                block_start = i - 1 # Convert to 0-based for UCSC browser
                block_height = height
            elif block_height == height:
                # Still part of the existing block
                continue
            else:
                # Convert to 0-based for UCSC, so this gives the interval (block_start, block_end)
                # which is half-open 0 based and hence goes up to but not including the 1-based position i
                block_end = i - 1

                # Don't output anything for 0s, but for everything else output the block
                if block_height != 0:
                    # Write out the block to the bed file
                    cov_file.write('\t'.join([chrom,
                                            str(block_start),
                                            str(block_end),
                                            str(block_height)]) + "\n")

                # Start a new block
                block_start = i - 1# Convert to 0-based for UCSC browser
                block_height = height

        if block_height is not None and block_height != 0:
            # Output the last block
                block_end = i - 1 # Convert to 0-based for UCSC browser
                # Write out the block to the bed file
                cov_file.write('\t'.join([chrom,
                                        str(block_start),
                                        str(block_end),
                                        str(block_height)]) + "\n")
