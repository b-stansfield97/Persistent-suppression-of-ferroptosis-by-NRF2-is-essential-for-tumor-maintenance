###############################################
##### Cell Marker Analysis and Assingment #####
###############################################

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

### Load in data from Data_Processing.r 
load("scRNA-Seq_initial_data_processing.rda")

############################
### Cell marker analysis ###
############################

### Convert data slot to RNA to run cluster marker analysis 
DefaultAssay(integrated_data) <- "RNA" 

### Normalize and Scale Raw Counts Data
integrated_data_RNA <- NormalizeData(integrated_data, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(integrated_data_RNA)
integrated_data_RNA <- ScaleData(integrated_data_RNA, features = all.genes)

### find  markers for every cell cluster in each sample
all_cluster_markers <- list()
for(i in unique(integrated_data_RNA@meta.data$sample)){
    data <- subset(integrated_data_RNA, subset = sample == i)
    cluster_markers2 <- FindAllMarkers(data, slot = "data", normalization.method = "LogNormalize", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    cluster_markers2$sample <- i
    all_cluster_markers[[i]] <- cluster_markers2
print(i)
}

### Collapse markers into one data frame
all_cluster_markers_comb <- data.table::rbindlist(all_cluster_markers)
all_cluster_markers_comb <- as.data.frame(all_cluster_markers_comb)

### Identify top 20 marker genes in each cluster 
top_20_cluster_markers <- list()
for(i in names(all_cluster_markers)){
    df <- all_cluster_markers[[i]]
    df <- df %>%
            group_by(cluster) %>%
            slice_max(n = 20, order_by = avg_log2FC)
    df$sample <- i
    top_20_cluster_markers[[i]] <- df
}

### Combine top 20 marker genes into one data frame and write data to CSV
top_20_cluster_markers_comb <- data.table::rbindlist(top_20_cluster_markers)
top_20_cluster_markers_comb <- as.data.frame(top_20_cluster_markers_comb)

write.csv(all_cluster_markers_comb, file = "all_cluster_markers.csv")
write.csv(top_20_cluster_markers_comb, file = "top_20_cluster_markers.csv")


### Plot cell markers to identify cell types
df <- read.csv("data/cell_markers.csv")
cell_markers <- split(df$marker, df$cell_type)

### Dotplots
for(i in names(cell_markers)){
    p1 <- DotPlot(subset(integrated_data_RNA, sample == "counts_ctrl"), features = cell_markers[[i]], cols = c("lightgrey", "blue"), , cluster.idents = TRUE)
    ggsave(p1, file = paste(i, "cell_markers.pdf", sep = "_"), height = 10, width = 10)
print(i)
}

### UMAPS
### Fuction to output UMAPs for vizual analysis of cell markers
plot_features <- function(genes, columns, figwidth, figheight, name){
    p2 <- FeaturePlot(integrated_data_RNA, features = genes, min.cutoff = "q9", split.by = "sample", ncol = columns, order = TRUE, slot = "scale.data", label = TRUE, cols = c("lightgrey", "red"))
    ggsave(p2, file = paste(name, "markers.pdf", sep = "_"), width = figwidth, height = figheight)
return(p2)
}

### PLot UMAPS
for (i in names(cell_markers)) {
    plot_features(genes     = cell_markers[[i]],
                  columns   = 1,
                  figwidth  = 16,
                  figheight = length(cell_markers[[i]]) * 4,
                  name      = i)
}


### Assign Cell Identities
# Cell type assignments are based on marker gene expression from the analysis above,
# cross-referenced with canonical markers from the literature for mouse lung cell types.

cell_types <- c("0" = "Endothelial",
                "1" = "Alveolar Type 2",
                "2" = "Alveolar Fibroblasts",
                "3" = "Mesenchymal",
                "4" = "Endothelial",
                "5" = "Aerocyte",
                "6" = "Epithelial",
                "7" = "Alveolar Type 1",
                "8" = "Endothelial",
                "9" = "Pericytes",
                "10" = "Epithelial",
                "11" = "Endothelial",
                "12" = "Endothelial",
                "13" = "Epithelial",
                "14" = "Club Cells",
                "15" = "Aerocytes",
                "16" = "Mesenchymal",
                "17" = "Lympatic",
                "18" = "Adventital Fibroblasts",
                "19" = "Alveolar Fibroblasts",
                "20" = "Pericytes",
                "21" = "Club Cells",
                "22" = "Mesenchymal",
                "23" = "Smooth Muscle",
                "24" = "Alveolar Type 2",
                "25" = "Alveolar Type 2",
                "26" = "Mesothelial",
                "27" = "Endothelial",
                "28" = "Epithelial")

### Assign cell types
integrated_data <- RenameIdents(integrated_data, cell_types)

### Save data
save(file = "Cell_Assigned_intergrated_data.rda")