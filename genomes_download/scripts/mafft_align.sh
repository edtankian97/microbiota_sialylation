#! /bin/bash

#Sequences alignment with mafft of sequence database

OUT_PATH=../Protein_database/mafft_align/
IN_PATH=../Protein_database/
mkdir -p $OUT_PATH

for file in $IN_PATH/*fasta; do
BASENAME=$(basename $file .fasta)

mafft --auto  --thread 5 $file > "$OUT_PATH/${BASENAME}_mafft.fasta"
done


