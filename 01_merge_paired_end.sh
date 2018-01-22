#!/bin/bash

# Add pear to path
export PATH=$PATH:/media/harddrive/tools/pear-0.9.11-linux-x86_64/bin/

# Create output dir
[ -d DNA ] || mkdir DNA

# run pear
pear -f data/DBC_S1_L001_R1_001.fastq.gz -r data/DBC_S1_L001_R2_001.fastq.gz -o DNA/merged.fastq -j 8

# compress everything
gzip DNA/*fastq
