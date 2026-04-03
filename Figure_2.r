###########################
######   Figure 2    ######
###########################

### Description
# R script to generate Figure 2 of the manuscript.
# Produces UMAP visualisations of scRNA-seq data following cell type assignment
# and box plots of cell cluster proportions for clusters of interest.

### Author
# Ben Stansfield

### Load in data after cell marker assignment
load("Cell_Assigned_intergrated_data.rda")

### Manually Plot UMAPs for better figure control ###
### get UMAP coordinates
UMAP <- as.data.frame(integrated_data[["umap"]]@cell.embeddings)

### get sample and cluster information for cells 
UMAP$sample <- integrated_data_RNA@meta.data$sample
UMAP$cluster <- integrated_data_RNA@meta.data$seurat_clusters
UMAP$clasifier <- integrated_data_RNA@meta.data$clasifier
UMAP$cellType <- Idents(integrated_data)

### asign levels for plotting order
UMAP$sample <- factor(UMAP$sample, levels = c("counts_ctrl","counts_2_1","counts_3_2","counts_4_3"))

### Plot sample specific UMAPs with colors assigned to cell types
### This plot was used for Figure 2B in the Manuscript

### Assign colors 
my_colors <- c("red", "blue", "darkgreen", 
                "orange", "purple", "#00FFFF", 
                "pink", "#FFFF00", "#A52A2A", 
                "tomato1", "black", "lightslateblue", 
                "#8B0000", "#808080")


annotated_umap <- ggplot(UMAP, aes(x = UMAP_1, y = UMAP_2, color = cellType, size = 0.01, alpha = 0.25))+
geom_point()+
scale_color_manual(values = my_colors)+
theme_classic(base_size = 12)+
scale_size_identity() + 
scale_alpha_identity() + 
guides(size = guide_legend(override.aes = list(size = 5)))+
facet_wrap(~sample, ncol = 4, nrow = 1)+
theme(strip.background = element_blank())
ggsave(annotated_umap, file = "Annoted_UMAP_by_sample.pdf", height = 2, width = 7)


########################################
### Cell cluster proportion analysis ###
########################################

# Based on visual inspection of the UMAPs, certain clusters seemingly decrease,
# with increased time from NRF2 deletion. Therefore we investigated the cell populations, 
# for these clusters across time. The figure produced is Figure 2D in the manuscript.

### Create a summary table of cluster counts per sample
cluster_sample_counts <- as.data.frame(table(

  integrated_data$seurat_clusters, 
  integrated_data$sample
))

### Rename the columns
colnames(cluster_sample_counts) <- c("Cluster", "Sample", "Count")

### Calculate the proportion of cells per cluster in each sample
cluster_sample_props <- as.data.frame(cluster_sample_counts) %>%
  group_by(Sample) %>%
  mutate(Proportion = Count / sum(Count) * 100)

cluster_sample_props <- as.data.frame(cluster_sample_props)

### Seperate out clusters of interest
clusters <- c("10","6","7","1","25","21","28")

bar_plotting <- cluster_sample_props[which(cluster_sample_props$Cluster %in% clusters), ]

# Create a stacked bar plot
bar_plotting$Sample <- factor(bar_plotting$Sample, levels = c("counts_ctrl", "counts_2_1", "counts_3_2", "counts_4_3"))

proportion_plot <- ggplot(bar_plotting, aes(x = Cluster, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(x = "Clusters", y = "Fraction of cluster per sample (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~Sample, ncol = 4, nrow = 1)
ggsave(proportion_plot, file = "proportion_bar_plot.pdf", width = 13, height = 3)


### Save Data
save(file = "Figure_2_data.rda")