#! /bin/bash

# Loop for each file with faa extension
 for arquivo in ../control_proteins/*.faa; do
    if [ -f "$arquivo" ]; then
        nome=$(basename "$arquivo")
        # Add name of file for each fasta header 
        sed -i "s/^>/>${nome} /" "$arquivo"
    fi
done
