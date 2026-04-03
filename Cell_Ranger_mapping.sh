#!/bin/bash
#
# Cell Ranger Read Alignment
# Aligns FASTQ files from four samples to the mm10 mouse reference genome in parallel.
# Reference: GRCm38/mm10 (refdata-gex-mm10-2020-A), CellRanger v7.1.0
#
# Author: Ben Stansfield

FASTQ_DIR="/data/8-19-23-BenCellRanger/FASTQ"
TRANSCRIPTOME="/data/refdata-gex-mm10-2020-A"

cellranger count --id=counts_1-Ctrl \
    --fastqs=${FASTQ_DIR} \
    --sample=1-Ctrl \
    --expect-cells=10000 \
    --transcriptome=${TRANSCRIPTOME} \
    --localcores=20 &

cellranger count --id=counts_2-1 \
    --fastqs=${FASTQ_DIR} \
    --sample=2-1_wk \
    --expect-cells=10000 \
    --transcriptome=${TRANSCRIPTOME} \
    --localcores=20 &

cellranger count --id=counts_3-2 \
    --fastqs=${FASTQ_DIR} \
    --sample=3-2_wk \
    --expect-cells=10000 \
    --transcriptome=${TRANSCRIPTOME} \
    --localcores=20 &

cellranger count --id=counts_4-3 \
    --fastqs=${FASTQ_DIR} \
    --sample=4-3_wk \
    --expect-cells=10000 \
    --transcriptome=${TRANSCRIPTOME} \
    --localcores=20 &

wait
echo "All CellRanger jobs completed."
