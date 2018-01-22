#!/bin/bash

for f in DNA/ID*.dna
do

# Get filenames
filename=$(basename "$f")
filename=${filename%.*}

echo " "
echo "$filename"
echo " "

python DNA-master/dna/dna_update.py -d DNA/$filename.dna

done
