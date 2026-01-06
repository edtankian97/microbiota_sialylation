#! /bin/bash

#Sequences alignment with mafft of sequence database

mkdir -p mafft_align

for file in ./*fasta; do
BASENAME=$(basename $file .fasta)

mafft --auto  --thread 5 $file > "${BASENAME}_mafft.fasta" 
mv "${BASENAME}_mafft.fasta" ./mafft_align
done


