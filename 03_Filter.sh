#!/bin/bash

file="DNA/merged.fastq.assembled.fastq.gz"

# Pullseq was downloaded from https://github.com/bcthomas/pullseq

# Filter out reads shorter than 117 nt and convert to fasta
pullseq -i $file --min 117 --convert > ${file}_filtered_temp

# Filter out reads longer than 117 nt (impossible to set min and max at same number)
pullseq -i ${file}_filtered_temp --max 117 -l 200 > ${file}_filtered.fasta

rm ${file}_filtered_temp

# Identify reads starting with A and ending with C
grep "^A" -B 1 --no-group-separator ${file}_filtered.fasta | grep "C$" --no-group-separator -B 1 > DNA/startA_endC.fasta

# Identify reads starting with A and ending with G
grep "^A" -B 1 --no-group-separator ${file}_filtered.fasta | grep "G$" --no-group-separator -B 1 > DNA/startA_endG.fasta

# Identify reads starting with T and ending with G
grep "^T" -B 1 --no-group-separator ${file}_filtered.fasta | grep "G$" --no-group-separator -B 1 > DNA/startT_endG.fasta

# Identify reads starting with T and ending with C
grep "^T" -B 1 --no-group-separator ${file}_filtered.fasta | grep "C$" --no-group-separator -B 1 > DNA/startT_endC.fasta

# Identify reads starting with G and ending with T
grep "^G" -B 1 --no-group-separator ${file}_filtered.fasta | grep "T$" --no-group-separator -B 1 > DNA/startG_endT.fasta

# Identify reads starting with G and ending with A
grep "^G" -B 1 --no-group-separator ${file}_filtered.fasta | grep "A$" --no-group-separator -B 1 > DNA/startG_endA.fasta

# Identify reads starting with C and ending with A
grep "^C" -B 1 --no-group-separator ${file}_filtered.fasta | grep "A$" --no-group-separator -B 1 > DNA/startC_endA.fasta

# Identify reads starting with C and ending with T
grep "^C" -B 1 --no-group-separator ${file}_filtered.fasta | grep "T$" --no-group-separator -B 1 > DNA/startC_endT.fasta

# Take reverse complement of the reads that start with G|T and end with T|G
fastarevcomp DNA/startG_endT.fasta > DNA/startG_endT_revcomp.fasta
fastarevcomp DNA/startG_endA.fasta > DNA/startG_endA_revcomp.fasta
fastarevcomp DNA/startC_endA.fasta > DNA/startC_endA_revcomp.fasta
fastarevcomp DNA/startC_endT.fasta > DNA/startC_endT_revcomp.fasta

# Remove the original files which are reverse complemented
rm -f DNA/startG_endT.fasta DNA/startG_endA.fasta DNA/startC_endA.fasta DNA/startC_endT.fasta

# Merge all files together and compress
cat DNA/start*fasta > ${file}_filtered.fasta_temp
rm -f DNA/start*

# Convert from multiline fasta to single line
fasta_formatter -i ${file}_filtered.fasta_temp -o ${file}_filtered.fasta -w 0
rm ${file}_filtered.fasta_temp

# Trim 1 nt in the beginning and 1 nt in the end to remove the orientation identifiers
fastx_trimmer -i ${file}_filtered.fasta -o ${file}_filtered_trimmed.fasta -f 2 -l 116

# Dereplicate
fastx_collapser -i ${file}_filtered_trimmed.fasta -o ${file}_filtered_trimmed_derep.fasta

# Convert to tsv
fasta_formatter -i ${file}_filtered_trimmed_derep.fasta -t > ${file}_filtered_trimmed_derep.tsv

# Compress filtered file
bzip2 ${file}_filtered.fasta
bzip2 ${file}_filtered_trimmed.fasta
bzip2 ${file}_filtered_trimmed_derep.fasta


