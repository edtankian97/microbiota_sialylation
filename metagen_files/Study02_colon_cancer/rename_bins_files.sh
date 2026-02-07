#!/bin/bash

for file in ./bins_paired/ERR*/bin*.fa
do
directory_name=$(dirname $file)
accession=$(basename $directory_name)
mv "${file}" "${directory_name}/${accession}_$(basename $file)"
done
