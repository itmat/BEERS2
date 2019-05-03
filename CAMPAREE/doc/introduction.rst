Introduction
============

CAMPAREE is a RNA expression simulator that is primed used real data to give realistic output.
CAMPREE needs as input a reference genome with transcript annotations as well as fastq files of samples of the species to base the output on.
For each sample, CAMPAREE outputs a simulated set of RNA transcripts mimicking expression levels with in the fastq files and accounting for isoform-level expression and allele-specific expression.
It also outputs simulated diploid genomes and their corresponding annotations with phased SNP and indel calls in the transcriptome from fastq reads.
Additionally the simulation outputs the underlying distributions used for expressing the transcripts.

The output transcripts are full-length RNA transcripts or pre-mRNA transcripts with introns included.


About
-----

CAMPAREE is prodcued at the Institute for Translational Medicine and Therapeutics, University of Pennsylvania.
