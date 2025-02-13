cellranger count --id=counts_1-Ctrl \
--fastqs=/data/8-19-23-BenCellRanger/FASTQ \
--sample=1-Ctrl \
--expect-cells=10000 \
--transcriptome=/data/refdata-gex-mm10-2020-A \
--localcores=20 &

cellranger count --id=counts_2-1 \
--fastqs=/data/8-19-23-BenCellRanger/FASTQ \
--sample=2-1_wk \
--expect-cells=10000 \
--transcriptome=/data/refdata-gex-mm10-2020-A \
--localcores=20 &

cellranger count --id=counts_3-2 \
--fastqs=/data/8-19-23-BenCellRanger/FASTQ \
--sample=3-2_wk \
--expect-cells=10000 \
--transcriptome=/data/refdata-gex-mm10-2020-A \
--localcores=20 &

cellranger count --id=counts_4-3 \
--fastqs=/data/8-19-23-BenCellRanger/FASTQ \
--sample=4-3_wk \
--expect-cells=10000 \
--transcriptome=/data/refdata-gex-mm10-2020-A \
--localcores=20 &
