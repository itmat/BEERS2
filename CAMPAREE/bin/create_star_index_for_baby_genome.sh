#!/bin/bash

# This bash script uses the STAR executable packaged with CAMPAREE to create
# a STAR index for the baby genome. This index is created in the
# resources/baby_genome.mm10/star_index.genome/ directory in the CAMPAREE
# path.

# Create directory for babY_gneome STAR index (if it doesn't already exist)
if [ ! -d "resources/baby_genome.mm10/star_index.genome" ]
then
    mkdir resources/baby_genome.mm10/star_index.genome
fi

#Run STAR to create the index (if it's destination directory was created correctly)
if [ -d "resources/baby_genome.mm10/star_index.genome" ]
then
third_party_software/STAR \
    --runMode genomeGenerate \
    --runThreadN 1 \
    --genomeDir resources/baby_genome.mm10/star_index.genome \
    --outFileNamePrefix resources/baby_genome.mm10/star_index.genome/ \
    --genomeFastaFiles resources/baby_genome.mm10/baby_genome.mm10.oneline_seqs.fa
else
    echo "There was a problem creating the directory for the STAR index."
fi
