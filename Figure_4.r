###########################
######   Figure 4    ######
###########################

### Description
# R script to generate Figure 4 of the manuscript.
# Produces:
# 1. Heatmaps of ferroptosis marker gene expression across cell clusters and time points.
# 2. Box plots of ferroptosis marker expression split by genotype at each time point.
# 3. UMAPs of ferroptosis marker expression across time points.

# Markers were selected based on experimental evidence that cells are undergoing
# ferroptosis following NRF2 deletion.

### Author
# Ben Stansfield

### Load data
load("post_classifier_analysis.rda")

### Plot heatmaps of ferroptosis marker genes - Figure 4A
ferroptosis_features <- c("Ptgs2", "Acsl4", "Tfrc", "Gapdh")

# Named vector mapping sample IDs to output file labels
samples <- c("counts_ctrl" = "ctrl",
             "counts_2_1"  = "week1",
             "counts_3_2"  = "week2",
             "counts_4_3"  = "week3")

for (i in names(samples)) {
  p <- DoHeatmap(subset(integrated_data_RNA, sample == i),
      features     = ferroptosis_features,
      size         = 3,
      slot         = "data",
      assay        = "RNA",
      group.colors = c("blue", "red", "green", "black"),
      group.by     = "clasifier",
      raster       = FALSE) +
    scale_fill_gradientn(colours  = c("white", "red", "black"),
                         limits   = c(0, 2),
                         na.value = "white")

  filename <- paste0(samples[[i]], "_NRF2_ferop_heatmap.pdf")
  ggsave(p, file = filename, height = 1, width = 10)
}


#################################################
### Ferroptosis Marker Proportion - Figure 4A ###
#################################################

# For each sample and each ferroptosis marker (Ptgs2, Acsl4, Tfrc), calculates
# the percentage of cells expressing the marker (count > 0) within each
# genotype class (classifier 0-3). This gives a per-genotype measure of how
# many cells are actively expressing each ferroptosis marker at each time point,
# allowing comparison of marker expression between NRF2 WT and KO cells over time.

### Fetch marker expression and genotype metadata
cell_express <- FetchData(integrated_data_RNA, vars = c("Ptgs2", "Acsl4", "Tfrc", "sample", "clasifier"))

cell_percentage <- list()
for(i in unique(cell_express$sample)){
    data <- cell_express[which(cell_express$sample == i), ]

    # get cell number for each genotype
    genotype_number <- table(data$clasifier)

    #Ptgs2
    Ptgs2_count <- table(data$Ptgs2, data$clasifier)
    Ptgs2_count_0_rem <- Ptgs2_count[which(rownames(Ptgs2_count) != "0"),]
    Ptgs2_ex_cel <- colSums(Ptgs2_count_0_rem)
    Ptgs2_percent <- (Ptgs2_ex_cel/genotype_number)*100
    names(Ptgs2_percent) <- paste("Ptgs2", i, names(Ptgs2_percent), sep = "-")
    cell_percentage[[i]] <- Ptgs2_percent

    # Acsl4
    Acsl4_count <- table(data$Acsl4, data$clasifier)
    Acsl4_count_0_rem <- Acsl4_count[which(rownames(Acsl4_count) != "0"),]
    Acsl4_ex_cel <- colSums(Acsl4_count_0_rem)
    Acsl4_percent <- (Acsl4_ex_cel/genotype_number)*100
    names(Acsl4_percent) <- paste("Acsl4", i, names(Acsl4_percent), sep = "-")
    cell_percentage[[i]] <- c(cell_percentage[[i]], Acsl4_percent)

    #Tfrc
    Tfrc_count <- table(data$Tfrc, data$clasifier)
    Tfrc_count_0_rem <- Tfrc_count[which(rownames(Tfrc_count) != "0"),]
    Tfrc_ex_cel <- colSums(Tfrc_count_0_rem)
    Tfrc_percent <- (Tfrc_ex_cel/genotype_number)*100
    names(Tfrc_percent) <- paste("Tfrc", i, names(Tfrc_percent), sep = "-")
    cell_percentage[[i]] <- c(cell_percentage[[i]], Tfrc_percent)
print(i)
}

## prep data for plotting
cell_percentage_plotting <- list()
for(i in names(cell_percentage)){
    data <- as.data.frame(cell_percentage[[i]])
    data$gene <- unlist(lapply(strsplit(rownames(data), "-"), function(x) x[1]))
    data$sample <- unlist(lapply(strsplit(rownames(data), "-"), function(x) x[2]))
    data$genotype <- unlist(lapply(strsplit(rownames(data), "-"), function(x) x[3]))
    colnames(data)[1] <- "percentage_of_cells"
    cell_percentage_plotting[[i]] <- data
}

### Combine all samples and plot bar charts of marker expression per genotype
library(data.table)
comb_data <- as.data.frame(rbindlist(cell_percentage_plotting))
comb_data$sample <- factor(comb_data$sample, levels = c("counts_ctrl", "counts_2_1", "counts_3_2", "counts_4_3"))

genotype_colors <- c("0" = "blue", "1" = "red", "2" = "green", "3" = "black")

for (gene in ferroptosis_features) {
  p <- ggplot(comb_data[which(comb_data$gene == gene), ],
              aes(x = sample, y = percentage_of_cells,
                  group = genotype, fill = genotype)) +
    geom_bar(stat = "identity", position = position_dodge(),
             color = "black", lwd = 0.5) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    scale_fill_manual(values = genotype_colors) +
    theme(legend.position = "none")

  ggsave(p, file = paste0(gene, "_cell_percentage_test.pdf"),
         height = 3, width = 6)
}


##################################################
### Ferroptosis Marker UMAPs - Figure 4B       ###
##################################################

# UMAPs of the integrated scRNA-seq data split by time point (control, 1-week,
# 3-week, 4-week post-NRF2 deletion), with cells coloured by expression level
# of each ferroptosis marker (Ptgs2, Acsl4, Tfrc). Cells not expressing the
# marker are shown in grey; expressing cells are coloured by expression level.
# This visualises the spatial distribution of ferroptosis marker expression
# across cell clusters and how it changes following NRF2 deletion.

### Plot ferroptosis marker UMAPs - 3 genes x 4 time points = 12 panels
integrated_data$sample <- factor(integrated_data$sample,
                                 levels = c("counts_ctrl", "counts_2_1",
                                            "counts_3_2", "counts_4_3"))

# Gapdh excluded here as it is a housekeeping control, not a ferroptosis marker
umap_genes <- c("Ptgs2", "Acsl4", "Tfrc")

ferop_plots <- list()
for (gene in umap_genes) {
  ferop_plots[[gene]] <- FeaturePlot(integrated_data,
                                     features  = gene,
                                     order     = TRUE,
                                     split.by  = "sample",
                                     pt.size   = 0.05,
                                     cols      = c("lightgrey", "red"))
}

# Combine into a single figure: 3 rows (genes) x 4 columns (time points)
ferop_combined <- ggarrange(plotlist = ferop_plots, nrow = 3, ncol = 1)
ggsave(ferop_combined, file = "ferroptosis_marker_UMAPs.pdf", height = 9, width = 12)

### Save data
save(file = "Final_data.rda")