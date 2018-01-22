#!/bin/bash

[ -d 02_QC_report_filtered ] || mkdir 02b_QC_report_filtered

fastqc DNA/merged.fastq.assembled.fastq.gz -o 02_QC_report
