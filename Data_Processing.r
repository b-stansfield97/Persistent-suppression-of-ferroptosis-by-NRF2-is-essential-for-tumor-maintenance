###################################
##### Initial Data processing #####
###################################

### Description
# R Script to read in counts, feautre and barcode data from scRNA-Seq experiment,
# assign cell genotypes, scTransform normalize, intergrate and perform dimensionality reduction 

### Author 
# Ben Stansfield

## Set working environment and load packages
rm(list=ls())
setwd("/data/8-19-23-BenCellRanger/FASTQ/Seurat")
library(Seurat)
library(ggplot2)
library(ggpubr)
Sys.setenv("R_MAX_VSIZE" = 1000000000)

## Pull file paths from directory
dirs <- list.dirs(path = "/data/8-19-23-BenCellRanger/FASTQ/Seurat/", recursive = F, full.names = F)
dirs <- dir(pattern = "counts")

# Generate seurat object for each sample 
seurat_list <- list()
for(i in dirs){
  name <- gsub("_filtered_feature_bc_matrix", "", i)
 cts <- ReadMtx(mtx = paste0("/data/8-19-23-BenCellRanger/FASTQ/Seurat/",i,"/matrix.mtx"),
                features = paste0("/data/8-19-23-BenCellRanger/FASTQ/Seurat/",i,"/features.tsv"),
                cells = paste0("/data/8-19-23-BenCellRanger/FASTQ/Seurat/",i,"/barcodes.tsv"))
seurat_list[[name]] <- assign(name, CreateSeuratObject(counts = cts))    
print(i)
}

### add sample name to seurat object meta data
for(i in names(seurat_list)){
    seurat_list[[i]]$sample <- rep(i, length.out = length(rownames(seurat_list[[i]]@meta.data)))
}

#####################################################
#### assign machine learning clasifiers to cells ####
#####################################################

# read in clasifier information from ML models
classifer <- list()
class <- dir(pattern = "*.txt") 
for(i in class){
    data <- read.delim2(i, header = FALSE)
    name <- gsub("*.txt", "", i)
    classifer[[name]] <- data
print(i)
}

## add clasifier information to meta data 
seurat_list[["counts_ctrl"]]$clasifier <- classifer[["counts_ctrl"]]$V1
seurat_list[["counts_2_1"]]$clasifier <- classifer[["counts_2-1"]]$V1
seurat_list[["counts_3_2"]]$clasifier <- classifer[["counts_3-2"]]$V1
seurat_list[["counts_4_3"]]$clasifier <- classifer[["counts_4-3"]]$V1

### remove cells which were classified differently by the two seperate models - clasifier 4
for(i in names(seurat_list)){
    sub_data  <- subset(seurat_list[[i]], clasifier != 4)
    seurat_list[[i]] <- sub_data
print(i)
}

###############################################################
#### QC - Data normalization and Scaling using SCTransform ####
###############################################################

### Normalization is run seperatly on each sample seperatly before intergrating the data

### get mitochondrial reads for each cell
for(i in names(seurat_list)){
    seurat_list[[i]]$percent.mt <- PercentageFeatureSet(seurat_list[[i]], pattern = '^mt-')
print(i)
}

### Normalize, Scale, find variable features and filter for percentage of mitochondrial genes
for(i in names(seurat_list)){
    seurat_list[[i]] <- SCTransform(seurat_list[[i]], vars.to.regress = "percent.mt", verbose = TRUE, variable.features.n = 5000, ncells = 10000)
print(i)
}

### 1-week post-Nrf2 deletion sample was used to determine anchor features because it had the highest representation of all four genotypes.
features <- FindVariableFeatures(seurat_list[["counts_2_1"]], selection.method = "vst", nfeatures = 5000)
features <- VariableFeatures(features)

### Variable features from 1-week post-Nrf2 deletion also present in every other sample following SCTransform were selected
common_features <- Reduce(intersect, lapply(seurat_list, function(x) rownames(x[["SCT"]]@scale.data)))
common_features_V2 <- intersect(common_features, features)

### assign common variable features to every sample
for(i in names(seurat_list)){
    VariableFeatures(seurat_list[[i]]) <- common_features_V2
print(i)
}

### run SCTintergration
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = common_features_V2)
cell_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = common_features_V2)
integrated_data <- IntegrateData(anchorset = cell_anchors, normalization.method = "SCT", verbose = TRUE)
DefaultAssay(integrated_data) <- "integrated"

## Dimension Reductionality
integrated_data <- RunPCA(integrated_data, verbose = TRUE, npcs = 50)
integrated_data <- FindNeighbors(integrated_data, reduction = "pca", dims = 1:20)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)
integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:20)
