# Persistent Suppression of Ferroptosis by NRF2 is Essential for Tumor Maintenance

scRNA-seq analysis pipeline for investigating NRF2-mediated ferroptosis suppression in a conditional Kras-driven mouse lung tumor model across multiple timepoints following Nrf2 deletion.

## Overview

This repository contains the bioinformatics pipeline used to process, integrate, and analyze single-cell RNA sequencing (scRNA-seq) data from mouse lung tissue. The study examines the role of NRF2 in tumor maintenance by profiling four distinct cell genotypes across four experimental timepoints following conditional Nrf2 deletion.

## Repository Structure

```
├── Cell_Ranger_mapping.sh   # CellRanger alignment of FASTQ reads to mm10 reference
├── Data_Processing.r        # Seurat-based QC, normalization, integration, and dimensionality reduction
└── README.md
```

## Requirements

### Software
- [CellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) v7.1.0
- R (≥ 4.0)
  - [Seurat](https://satijalab.org/seurat/) v4+
  - ggplot2
  - ggpubr

### Reference Genome
- Mouse reference genome: GRCm38/mm10 (`refdata-gex-mm10-2020-A`)

## Pipeline

### 1. Read Alignment (`Cell_Ranger_mapping.sh`)
FASTQ files from four samples are aligned to the mm10 mouse reference genome using CellRanger. All four samples are processed in parallel:

| Sample ID   | Timepoint            |
|-------------|----------------------|
| counts_ctrl | Control (corn oil)   |
| counts_2_1  | 1 week post-Nrf2 KO  |
| counts_3_2  | 2 weeks post-Nrf2 KO |
| counts_4_3  | 3 weeks post-Nrf2 KO |

### 2. Data Processing (`Data_Processing.r`)
1. Load CellRanger count matrices into Seurat objects
2. Assign cell genotypes using pre-trained ML classifiers
3. QC filtering and SCTransform normalization per sample
4. Data integration using Seurat's `IntegrateData` (SCT workflow)
5. Dimensionality reduction: PCA → UMAP

## Experimental Design

To investigate the role of NRF2 in tumor progression within the tumor microenvironment, we generated a conditional Nrf2 knockout mouse model in a Kras<sup>G12D</sup>-driven tumor system: *Kras*<sup>FSF.G12D/+</sup>;*Nrf2*<sup>Fl/Fl</sup>;*Rosa26*<sup>CreERT2/CreERT2</sup> (KNR) mice. 8-week-old mice were instilled with FlpO virus to activate Kras<sup>G12D</sup> for tumor initiation. Mice were then treated with tamoxifen (TAM) to conditionally knockout Nrf2 for 1, 2, or 3 weeks prior to tissue harvest. Control mice were instilled with FlpO virus but received corn oil instead of TAM.

## Cell Genotyping

The model contains cells with four distinct genotypes (*Kras*<sup>+/+</sup>;*Nrf2*<sup>+/+</sup>, *Kras*<sup>+/+</sup>;*Nrf2*<sup>-/-</sup>, *Kras*<sup>G12D/+</sup>;*Nrf2*<sup>+/+</sup>, *Kras*<sup>G12D/+</sup>;*Nrf2*<sup>-/-</sup>). Cell genotypes were assigned using binary classification neural network models trained on cells with complete alignment to the Kras locus and reliable Nrf2 exon 5 coverage, allowing clear distinction between Kras<sup>+/+</sup> and Kras<sup>G12D/+</sup>, and between Nrf2<sup>+/+</sup> and Nrf2<sup>-/-</sup> cells.

## Data

scRNA-seq libraries were prepared using the 10x Genomics Chromium Single Cell 3' Reagent Kits v3.1 and sequenced on an Illumina NovaSeq S4 (2×100 bp). Demultiplexing was performed with Bcl2fastq (v2.19.1.403).

## License

This project is licensed under the terms of the [LICENSE](LICENSE) file included in this repository.
