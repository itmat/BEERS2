Introduction
============

BEERS2 is an advanced RNA-seq simulator that simulates invidual molecules throughout the entire library preparation and sequencing process.
It is designed to operate on the output of `CAMPAREE <https://github.com/itmat/CAMPAREE>`_, which simulates realistic samples of RNA for input to sequencing.

BEERS2 supports various biases present in real sequencing, including GC bias, primer sequencing bias, positional biases (from polyA selection), and RiboZero depletion bias, as well as typical read errors.
BEERS2 can output as FASTQ files or as SAM/BAM files with ideal alignments of the reads.
This makes it ideal for benchmarking of downstream analyses in the RNA-seq pipeline, such as quantification.
