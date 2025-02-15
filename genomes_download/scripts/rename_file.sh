#!/bin/bash

for file in ../proteins/ncbi_dataset/data/*/protein.faa
do
directory_name=$(dirname $file)
accession=$(basename $directory_name)
mv "${file}" "${directory_name}/${accession}_$(basename $file)"
done
