###################################
##### Initial Data Processing #####
###################################

### Description
# R script to read in counts, features and barcode data from scRNA-seq experiment,
# assign cell genotypes, SCTransform normalize, integrate and perform dimensionality reduction.

### Author
# Ben Stansfield

## Set working environment and load packages
rm(list = ls())
setwd("/data/8-19-23-BenCellRanger/FASTQ/Seurat")
library(Seurat)
library(ggplot2)
library(ggpubr)
Sys.setenv("R_MAX_VSIZE" = 1000000000)

## Pull file paths from directory
dirs <- list.dirs(path = "/data/8-19-23-BenCellRanger/FASTQ/Seurat/", recursive = F, full.names = F)
dirs <- dir(pattern = "counts")

# Generate Seurat object for each sample
seurat_list <- list()
for (i in dirs) {
    name <- gsub("_filtered_feature_bc_matrix", "", i)
    cts <- ReadMtx(mtx      = paste0("/data/8-19-23-BenCellRanger/FASTQ/Seurat/", i, "/matrix.mtx"),
                   features = paste0("/data/8-19-23-BenCellRanger/FASTQ/Seurat/", i, "/features.tsv"),
                   cells    = paste0("/data/8-19-23-BenCellRanger/FASTQ/Seurat/", i, "/barcodes.tsv"))
    seurat_list[[name]] <- assign(name, CreateSeuratObject(counts = cts))
    print(i)
}

### Add sample name to Seurat object metadata
for (i in names(seurat_list)) {
    seurat_list[[i]]$sample <- rep(i, length.out = length(rownames(seurat_list[[i]]@meta.data)))
}

######################################################
#### Assign Machine Learning Classifiers to Cells ####
######################################################

# Read in classifier information from ML models
classifier <- list()
classifier_files <- dir(pattern = "*.txt")
for (i in classifier_files) {
    data <- read.delim2(i, header = FALSE)
    name <- gsub("*.txt", "", i)
    classifier[[name]] <- data
    print(i)
}

## Add classifier information to metadata
seurat_list[["counts_ctrl"]]$classifier <- classifier[["counts_ctrl"]]$V1
seurat_list[["counts_2_1"]]$classifier  <- classifier[["counts_2-1"]]$V1
seurat_list[["counts_3_2"]]$classifier  <- classifier[["counts_3-2"]]$V1
seurat_list[["counts_4_3"]]$classifier  <- classifier[["counts_4-3"]]$V1

### Remove cells classified differently by the two separate models (classifier 4)
for (i in names(seurat_list)) {
    sub_data <- subset(seurat_list[[i]], classifier != 4)
    seurat_list[[i]] <- sub_data
    print(i)
}

###############################################################
#### QC - Data Normalization and Scaling using SCTransform ####
###############################################################

### Normalization is run separately on each sample before integrating the data

### Get mitochondrial reads for each cell
for (i in names(seurat_list)) {
    seurat_list[[i]]$percent.mt <- PercentageFeatureSet(seurat_list[[i]], pattern = '^mt-')
    print(i)
}

### Normalize, scale, find variable features and filter for percentage of mitochondrial genes
for (i in names(seurat_list)) {
    seurat_list[[i]] <- SCTransform(seurat_list[[i]], vars.to.regress = "percent.mt", verbose = TRUE, variable.features.n = 5000, ncells = 10000)
    print(i)
}

### The 1-week post-Nrf2 deletion sample was used to determine anchor features
### as it had the highest representation of all four genotypes.
features <- FindVariableFeatures(seurat_list[["counts_2_1"]], selection.method = "vst", nfeatures = 5000)
features <- VariableFeatures(features)

### Variable features from the 1-week post-Nrf2 deletion sample also present in
### every other sample following SCTransform were selected.
common_features    <- Reduce(intersect, lapply(seurat_list, function(x) rownames(x[["SCT"]]@scale.data)))
common_features_V2 <- intersect(common_features, features)

### Assign common variable features to every sample
for (i in names(seurat_list)) {
    VariableFeatures(seurat_list[[i]]) <- common_features_V2
    print(i)
}

### Run SCT integration
seurat_list     <- PrepSCTIntegration(object.list = seurat_list, anchor.features = common_features_V2)
cell_anchors    <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = common_features_V2)
integrated_data <- IntegrateData(anchorset = cell_anchors, normalization.method = "SCT", verbose = TRUE)
DefaultAssay(integrated_data) <- "integrated"

## Dimensionality Reduction
integrated_data <- RunPCA(integrated_data, verbose = TRUE, npcs = 50)
integrated_data <- FindNeighbors(integrated_data, reduction = "pca", dims = 1:20)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)
integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:20)

### Plot UMAP
DimPlot(integrated_data, reduction = "umap", split.by = "sample")

### Save Data
save.image(file = "scRNA-Seq_initial_data_processing.rda")
