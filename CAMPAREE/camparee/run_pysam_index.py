import sys
import os
import pysam

"""Simple command line wrapper around pysam.index()."""

if len(sys.argv) > 1:
    bam_filename = sys.argv[1]
    if os.path.isfile(bam_filename):
        pysam.index(bam_filename)
    else:
        print(f"ERROR: {bam_filename} cannot be accessed or does not exist.", file=sys.stderr)
else:
    print(f"ERROR: no bam filename provided.", file=sys.stderr)
