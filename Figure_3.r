###########################
######   Figure 3    ######
###########################

### Description
# R script to generate Figure 3 of the manuscript.
# Produces 
# 1. Visualizations of the example stratergy to determine cell genotypes.
# 2. Alluvial Diagram demostrating the changes in cell genotype over time 
# 3. UMAPs seperated by cell genotype over time

### Author
# Ben Stansfield

### Load libraries for BAM file visualisation
library(Rsamtools)
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(biomaRt)
library("BSgenome.Mmusculus.UCSC.mm10")
options(ucscChromosomeNames = FALSE)

###########################################
### NRF2 WT vs KO BAM File Preparation  ###
###########################################

# The ML classifier outputs cell barcodes without the Cell Ranger GEM group suffix.
# The "-1" suffix is appended here to match the barcode format in the BAM file (e.g. ACGTACGTACGT-1).
# The updated barcode lists are written to new text files for downstream BAM file generation.

### Append "-1" suffix to NRF2 WT barcodes
wt <- read.delim2("c2-1-nrf2wt.txt", header = FALSE)
wt$V1 <- paste(wt$V1, "-1", sep = "")
sink(file = "Updated_c2-1-nrf2wt.txt")
cat(unlist(wt), sep="
")
sink()

### Append "-1" suffix to NRF2 KO barcodes
ko <- read.delim2("c2-1-nrf2ko.txt", header = FALSE)
ko$V1 <- paste(ko$V1, "-1", sep = "")
sink(file = "Updated_c2-1-nrf2ko.txt")
cat(unlist(ko), sep="
")
sink()

### Samtools BAM File Generation of NRF2 WT and KO Cells
# The following commands were run in the terminal (bash) to subset the BAM files
# by cell barcode, isolating NRF2 WT and KO cells for visualisation.
#
# [bash]
#
#   # Generate BAM file for NRF2 KO cells
#   samtools view -b possorted_genome_bam.bam -D CB:Updated_c2-1-nrf2ko.txt chr2 >> week1_nrf2_ko_bam.bam
#   samtools index week1_nrf2_ko_bam.bam
#
#   # Generate fusion BAM file for NRF2 KO cells (filters for reads containing the deletion junction)
#   samtools view -b -h -D CB:Updated_c2-1-nrf2ko.txt -o week1_nrf2_ko_fusion.bam possorted_genome_bam.bam chr2:75673763-75679881 | awk '/TGAGGGCAGGCAGAT/'

### Path to BAM files
Control_nrf2_wt <- "/Users/benstansfield/Desktop/NRF2_scRNA-Seq/BAM Analysis/Control/possorted_genome_bam.bam"
week1_ko_nrf2   <- "/Users/benstansfield/Desktop/NRF2_scRNA-Seq/BAM Analysis/Week 1/week1_nrf2_ko_fusion.bam"

### Assign Gviz tracks for NRF2 locus (chr2)
gtrack        <- GenomeAxisTrack()
itrack        <- IdeogramTrack(genome = "mm10", chromosome = "chr2")
txdb          <- TxDb.Mmusculus.UCSC.mm10.knownGene
grtrack       <- GeneRegionTrack(txdb, genome = "mm10", chromosome = "chr2", name = "Nrf2", start = 75675513, end = 75677668, fill = "orange")
control_allign     <- AlignmentsTrack(Control_nrf2_wt, genome = "mm10", chromosome = "chr2", name = "Nrf2+/+")
week_1_mt_allign   <- AlignmentsTrack(week1_ko_nrf2,   genome = "mm10", chromosome = "chr2", name = "Nrf2-/-")

### Plot and save NRF2 WT and KO tracks
pdf(file = "Control_Nrf2_wt.pdf", width = 5, height = 5)
plotTracks(list(itrack, grtrack, control_allign), from = 75675513 - 1000, to = 75677668 + 1000)
dev.off()

pdf(file = "Week1_Nrf2_ko.pdf", width = 5, height = 15)
plotTracks(list(grtrack, week_1_mt_allign), from = 75675513 - 1000, to = 75677668 + 1000)
dev.off()

###########################################################
### KRAS WT vs MT BAM File Preparation Figure 3 A and B ###
###########################################################

# The same barcode suffix correction is applied to the KRAS classifier outputs.

### Append "-1" suffix to KRAS WT barcodes
kras_wt <- read.delim2("c-2-1-kraswt.txt", header = FALSE)
kras_wt$V1 <- paste(kras_wt$V1, "-1", sep = "")
sink(file = "c-2-1-kraswt.txt")
cat(unlist(kras_wt), sep="
")
sink()

### Append "-1" suffix to KRAS MT barcodes
kras_mt <- read.delim2("c-2-1-krasmt.txt", header = FALSE)
kras_mt$V1 <- paste(kras_mt$V1, "-1", sep = "")
sink(file = "c-2-1-krasmt.txt")
cat(unlist(kras_mt), sep="
")
sink()

### Samtools BAM File Generation of KRAS WT and MT Cells
# The following commands were run in the terminal (bash) to subset the BAM files
# by cell barcode, isolating KRAS WT and G12D mutant cells for visualisation.
#
# [bash]
#
#   samtools view -b possorted_genome_bam.bam -D CB:c-2-1-kraswt.txt chr6:145246760-145246781 >> week1_kras_wt_bam.bam
#   samtools view -b possorted_genome_bam.bam -D CB:c-2-1-krasmt.txt chr6:145246760-145246781 >> week1_kras_mt_bam.bam

### Path to BAM files
week1_wt_kras_bamfile <- "/Users/benstansfield/Desktop/NRF2_scRNA-Seq/BAM Analysis/Week 1/week1_kras_wt_bam.bam"
week1_mt_kras_bamfile <- "/Users/benstansfield/Desktop/NRF2_scRNA-Seq/BAM Analysis/Week 1/week1_kras_mt_bam.bam"

### Assign Gviz tracks for KRAS locus (chr6)
gtrack        <- GenomeAxisTrack()
itrack        <- IdeogramTrack(genome = "mm10", chromosome = "chr6")
txdb          <- TxDb.Mmusculus.UCSC.mm10.knownGene
seq_track     <- SequenceTrack(BSgenome.Mmusculus.UCSC.mm10, chromosome = "chr6", from = 145246760 - 3, to = 145246781 + 3, cex = 0.5, min.height = 50)
week1_wt_kras <- AlignmentsTrack(week1_wt_kras_bamfile, isPaired = TRUE, genome = "mm10", type = "pileup", chromosome = "chr6", name = "Kras+/+")
week1_mt_kras <- AlignmentsTrack(week1_mt_kras_bamfile, isPaired = TRUE, genome = "mm10", type = "pileup", chromosome = "chr6", name = "KrasG12D/+")

### Plot and save KRAS WT and MT tracks
pdf(file = "Week1_Nrf2_wt_Nrf2_Mt.pdf", width = 3, height = 4)
plotTracks(list(week1_wt_kras, week1_mt_kras, seq_track), from = 145246760 - 3, to = 145246781 + 3)
dev.off()

##################################
### Alluvial Plot - Figure 3C  ###
##################################

# Alluvial plot showing how the proportion of cells in each genotype class
# (NRF2 WT, NRF2 KO, KRAS WT, KRAS MT) changes across time points following
# NRF2 deletion. Each flow represents the shift in cell genotype composition
# between the control, 1-week, 3-week, and 4-week post-deletion samples.

# For each sample, a contingency table is built counting how many cells in each
# cluster (column) belong to each ML classifier class (row: 0 = NRF2 WT/KRAS WT,
# 1 = NRF2 KO/KRAS WT, 2 = NRF2 WT/KRAS MT, 3 = NRF2 KO/KRAS MT).
# Counts are then normalised to the total number of cells in that sample to give
# the percentage of each genotype class per cluster per sample.

### Build per-sample contingency tables of classifier proportions
contingency_tables_proportion <- list()

for (i in unique(integrated_data@meta.data$sample)) {

    # Initialise empty matrix: rows = classifier classes, columns = clusters
    table <- matrix(0,
                    ncol = length(unique(integrated_data@meta.data$seurat_clusters)),
                    nrow = length(unique(integrated_data@meta.data$sample)))
    colnames(table) <- unique(integrated_data@meta.data$seurat_clusters)
    rownames(table) <- c(paste(i, "0", sep = "-"),
                         paste(i, "1", sep = "-"),
                         paste(i, "2", sep = "-"),
                         paste(i, "3", sep = "-"))

    # Subset to current sample and count classifier classes per cluster
    data <- subset(integrated_data, subset = sample == i)
    for (y in unique(data@meta.data$seurat_clusters)) {
        data2 <- data@meta.data[which(data@meta.data$seurat_clusters == y), ]
        table[1, y] <- length(which(data2$clasifier == "0"))
        table[2, y] <- length(which(data2$clasifier == "1"))
        table[3, y] <- length(which(data2$clasifier == "2"))
        table[4, y] <- length(which(data2$clasifier == "3"))
    }

    # Normalise counts to percentage of total cells in the sample
    table2      <- apply(table, 2, as.numeric)
    cell_number <- sum(table2)
    table3      <- apply(table2, 2, function(x) (x / cell_number) * 100)

    contingency_tables_proportion[[i]] <- table3
    print(i)
}

### Add sample and classifier metadata columns to each table
for (i in names(contingency_tables_proportion)) {
    data2 <- as.data.frame(contingency_tables_proportion[[i]])
    data2$sample     <- i
    data2$classifier <- seq(0, 3, 1)
    contingency_tables_proportion[[i]] <- data2
}

### Combine all samples into one data frame
contingency_tables_prop_comb <- data.table::rbindlist(contingency_tables_proportion)

### Get cluster-free genotype counts for alluvial plot
classifier_only <- contingency_tables_prop_comb
classifier_only <- classifier_only[, -c("sample", "classifier")]
summed_classifier <- data.frame(count = rowSums(classifier_only),
                                sample = contingency_tables_prop_comb$sample,
                                classifier = contingency_tables_prop_comb$classifier)

### Alluvial plots ###
library(alluvial)
library(ggalluvial)

summed_classifier$classifier <- as.factor(summed_classifier$classifier)
summed_classifier$sample <- as.factor(summed_classifier$sample)
summed_classifier$sample<- factor(summed_classifier$sample, levels = c("counts_ctrl", "counts_2_1", "counts_3_2", "counts_4_3"))

colors <- c("blue", "red", "green", "black")

alluvial_p <- ggplot(summed_classifier,
       aes(x = sample, alluvium = classifier,
           y = count,
           fill = classifier, label = classifier)) +
  geom_flow(aes(fill = classifier, color = classifier), width = 0, alpha = 0.75)+
  scale_fill_manual(values = c("0" = "blue", "1" = "red", "2" = "green", "3" = "black")) +
  scale_color_manual(values = c("0" = "blue", "1" = "red", "2" = "green", "3" = "black")) +
  theme_classic()+
  ylab("% of Total Cells")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))
       
ggsave(alluvial_p, file = "Classifier_time_alluvial.pdf", width = 7, height = 6)


######################################
### UMAP by Classifier - Figure 3D ###
######################################

# 16 UMAPs arranged as a 4x4 grid: one panel per genotype class (NRF2 WT/KRAS WT,
# NRF2 KO/KRAS WT, NRF2 WT/KRAS MT, NRF2 KO/KRAS MT), each containing 4 plots
# split by time point (control, 1-week, 3-week, 4-week post-NRF2 deletion).
# Cells are coloured by cluster identity, allowing visualisation of how the
# cluster composition of each genotype class shifts across time.

### Manually Plot UMAPs for better figure control ###
### get UMAP coordinates
UMAP <- as.data.frame(integrated_data[["umap"]]@cell.embeddings)

### get sample and cluster information for cells 
UMAP$sample <- integrated_data_RNA@meta.data$sample
UMAP$cluster <- integrated_data_RNA@meta.data$seurat_clusters
UMAP$clasifier <- integrated_data_RNA@meta.data$clasifier
UMAP$cellType <- Idents(integrated_data)

### get coordinates to plot cluster number
label_coord <- matrix(0, ncol = 2, nrow = length(unique(UMAP$cluster)))
colnames(label_coord) <- c("x","y")
rownames(label_coord) <- unique(UMAP$cluster)

for(i in unique(UMAP$cluster)){
    df <- UMAP[which(UMAP$cluster == i), ]
    label_coord[i, "x"] <- median(df$UMAP_1) # get median x coordinate
    label_coord[i, "y"] <- median(df$UMAP_2) # get median y coordinate
}

label_coord <- as.data.frame(label_coord)
label_coord <- label_coord[as.character(seq(0,28,1)), ]
label_coord$cellType <- as.character(cell_types)
label_coord$cluster <- rownames(label_coord)

### asign levels for plotting order
UMAP$sample <- factor(UMAP$sample, levels = c("counts_ctrl","counts_2_1","counts_3_2","counts_4_3"))

### Assign colors 
my_colors <- c("red", "blue", "darkgreen", 
                "orange", "purple", "#00FFFF", 
                "pink", "#FFFF00", "#A52A2A", 
                "tomato1", "black", "lightslateblue", 
                "#8B0000", "#808080")

### Generate one UMAP panel per genotype class (classifier 0-3), faceted by time point
# classifier 0 = NRF2 WT / KRAS WT
# classifier 1 = NRF2 KO / KRAS WT
# classifier 2 = NRF2 WT / KRAS MT
# classifier 3 = NRF2 KO / KRAS MT

classifier_plots <- list()

for (i in 0:3) {
    data <- subset(UMAP, clasifier == as.character(i))
    classifier_plots[[i + 1]] <- ggplot(data, aes(x = UMAP_1, y = UMAP_2, color = cellType, size = 0.01, alpha = 0.5)) +
        geom_point() +
        scale_color_manual(values = my_colors) +
        theme_classic(base_size = 12) +
        scale_size_identity() +
        scale_alpha_identity() +
        guides(size  = guide_legend(override.aes = list(size = 5)),
               alpha = guide_legend(override.aes = list(alpha = 1))) +
        geom_text(label_coord, mapping = aes(x = x, y = y, label = cluster), size = 4, color = "black", alpha = 1) +
        facet_wrap(~sample, ncol = 4, nrow = 1) +
        theme(strip.background = element_blank(),
              legend.position  = "none")
}

p6 <- ggarrange(plotlist = classifier_plots, align = "hv", ncol = 1, nrow = 4)


### Save data
save(file = "post_classifier_analysis.rda")