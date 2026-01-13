#!/bin/bash

# Loop for each .faa file extension
for arquivo in ./ERR*/*prokka/*.faa; do
    # Verify if file exist 
    if [ -f "$arquivo" ]; then
        # Add filename to each fasta header 
        sed -i "s/^>/>$arquivo /g" "$arquivo"
    fi
done
