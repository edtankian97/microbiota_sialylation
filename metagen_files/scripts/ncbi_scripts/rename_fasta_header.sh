#!/bin/bash

# Loop for each .faa file extension
 for arquivo in ./all_prokka/*.faa; do
    if [ -f "$arquivo" ]; then
        name=$(basename "$arquivo")
        sed -i "s/^>/>${name} /" "$arquivo"
    fi
done
