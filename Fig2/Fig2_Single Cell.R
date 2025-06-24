rm(list = ls())
.libPaths(c("~/SeuratV5", .libPaths()))

#---------------------#
# Step 0: Load packages
#---------------------#
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(patchwork)
  library(cowplot)
  library(RColorBrewer)
  library(harmony)
  library(biomaRt)
})

#---------------------#
# Step 1: Load Data
#---------------------#
load_10X_data <- function(path, label) {
  counts <- Read10X(data.dir = path)
  CreateSeuratObject(counts = counts, project = label)
}

spleen_age_path <- "/home/data/t160407/Organ_aging/Mouse Aging/Spleen/Row data/GSM4321530_Spleen_Aged"
spleen_young_path <- "/home/data/t160407/Organ_aging/Mouse Aging/Spleen/Row data/GSM4321531_Spleen_Young"

age_obj <- load_10X_data(spleen_age_path, "Age")
young_obj <- load_10X_data(spleen_young_path, "Young")

# Merge datasets
sce.all <- merge(x = age_obj, y = young_obj, add.cell.ids = c("Age", "Young"), project = "Mouse_spleen")
sce.all <- JoinLayers(sce.all)

#---------------------#
# Step 2: QC filtering
#---------------------#
sce.all$log10GenesPerUMI <- log10(sce.all$nFeature_RNA) / log10(sce.all$nCount_RNA)
sce.all$mitoRatio <- PercentageFeatureSet(sce.all, pattern = "^mt-") / 100

# Filter cells
sce.filtered <- subset(sce.all,
                       subset = nCount_RNA >= 500 &
                         nFeature_RNA >= 250 & nFeature_RNA <= 3000 &
                         log10GenesPerUMI > 0.8 &
                         mitoRatio < 0.15)

#---------------------#
# Step 3: Normalization
#---------------------#
sce.filtered <- NormalizeData(sce.filtered)

#---------------------#
# Step 4: Cell cycle scoring (mouse genes)
#---------------------#
convert_to_mouse <- function(human_genes) {
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  getLDS(attributes = "hgnc_symbol", filters = "hgnc_symbol", values = human_genes,
         mart = human, attributesL = "mgi_symbol", martL = mouse, uniqueRows = TRUE)$MGI.symbol
}

s_genes_mouse <- convert_to_mouse(cc.genes$s.genes)
g2m_genes_mouse <- convert_to_mouse(cc.genes$g2m.genes)

sce.filtered <- CellCycleScoring(sce.filtered, s.features = s_genes_mouse, g2m.features = g2m_genes_mouse)

#---------------------#
# Step 5–6: HVGs, Scaling, PCA
#---------------------#
sce.filtered <- FindVariableFeatures(sce.filtered, nfeatures = 2000)
sce.filtered <- ScaleData(sce.filtered)
sce.filtered <- RunPCA(sce.filtered)

#---------------------#
# Step 7: Harmony integration
#---------------------#
sce.filtered <- RunHarmony(sce.filtered, group.by.vars = "orig.ident")
sce.filtered <- RunUMAP(sce.filtered, reduction = "harmony", dims = 1:30)
sce.filtered <- FindNeighbors(sce.filtered, reduction = "harmony", dims = 1:30)
sce.filtered <- FindClusters(sce.filtered, resolution = 0.5)

#---------------------#
# Step 8: Visualization
#---------------------#
ncluster <- length(unique(sce.filtered$seurat_clusters))
mycol <- colorRampPalette(brewer.pal(12, "Set3"))(ncluster)

# UMAP
pdf("Spleen/Result/umap_cluster.pdf", width = 10, height = 6)
plot1 <- DimPlot(sce.filtered, reduction = "umap", label = TRUE, cols = mycol)
plot2 <- DimPlot(sce.filtered, reduction = "umap", group.by = "orig.ident")
print(plot1 + plot2)
dev.off()

# TSNE
sce.filtered <- RunTSNE(sce.filtered, reduction = "harmony", dims = 1:30)
pdf("Spleen/Result/tsne_cluster.pdf", width = 10, height = 6)
tsne1 <- DimPlot(sce.filtered, reduction = "tsne", label = TRUE, cols = mycol)
tsne2 <- DimPlot(sce.filtered, reduction = "tsne", group.by = "orig.ident")
print(tsne1 + tsne2)
dev.off()

#..................................................

####Calculating gene scores for aging trends######################
# Load libraries
library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(Seurat)
library(tidydr)
library(ggthemes)
library(reshape2)
library(RColorBrewer)

# Normalize vector to [0,1]
normalize_vector <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Plot AUC on t-SNE
plot_AUC_tsne <- function(pbmc, output_path) {
  p <- FeaturePlot(pbmc, features = "AUC", reduction = "tsne") +
    scale_color_gradientn(values = quantile(pbmc$AUC),
                          colours = c('#68879D', 'white', '#BC3C29FF')) +
    tidydr::theme_dr() +
    theme(panel.grid = element_blank()) +
    coord_fixed()
  ggsave(plot = p, filename = output_path, height = 5, width = 7)
}

# Plot top N% high AUC cells
plot_top_cells_bar <- function(meta_data, output_path) {
  p <- ggplot(meta_data, aes(x = group, fill = group)) +
    geom_bar(position = "dodge") +
    scale_fill_manual(values = c("#d1d2d2", "#fbd3b9", "#a1c9e5", "#417bb9")) +
    ylab("Cell numbers") +
    xlab("Top 10% cells expressing trend genes") +
    theme_par()
  ggsave(output_path, plot = p, height = 4, width = 6)
}

# Plot aging cell proportion per category
plot_aging_ratio <- function(top_cells, full_meta, output_path) {
  aging_ratio <- table(top_cells$category) / table(full_meta$category)
  df <- data.frame(
    category = names(table(full_meta$category)),
    aging_ratio = as.numeric(aging_ratio),
    non_aging_ratio = 1 - as.numeric(aging_ratio)
  )
  df <- melt(df, id.vars = "category")
  df$variable <- factor(df$variable, levels = c("non_aging_ratio", "aging_ratio"))
  df$category <- factor(df$category, levels = names(sort(aging_ratio, decreasing = TRUE)))
  
  p <- ggplot(df, aes(x = category, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack", color = "black", width = 0.7, size = 0.25) +
    scale_fill_manual(values = c("#68879D", "#BC3C29FF")) +
    labs(x = "", y = "Cell number") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white", colour = "black", size = 0.25),
      axis.line = element_line(colour = "black", size = 0.25),
      axis.title = element_text(size = 13, color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      legend.position = c(0.9, 0.85)
    )
  ggsave(plot = p, filename = output_path, height = 5, width = 7)
}

# ==== Main Analysis ====
{
  pbmc <- sce.filtered
  
  # Load gene set and categorize
  gene_data <- readRDS("/home/data/t140311/Organ_aging/gene_set/cluster_spleen_mfuzz_file_gene_exp.rds")
  gene_data$up_down <- ifelse(gene_data$mk.statistic > 0, "Up", "Down")
  gene_data <- gene_data[, c("gene", "up_down")]
  
  # Build gene sets
  gene_sets <- split(gene_data$gene, gene_data$up_down)
  
  # AUCell scoring
  cell_rankings <- AUCell_buildRankings(pbmc@assays$RNA@data, splitByBlocks = TRUE)
  auc_result <- AUCell_calcAUC(gene_sets, cell_rankings, aucMaxRank = nrow(cell_rankings) * 0.1)
  
  auc_up <- as.numeric(getAUC(auc_result)["Up", ])
  auc_down <- as.numeric(getAUC(auc_result)["Down", ])
  pbmc$AUC <- normalize_vector(auc_up - auc_down)
  
  # Plot AUC map
  plot_AUC_tsne(pbmc, "AUC.pdf")
  
  # Top 10% cells
  meta_data <- pbmc@meta.data
  top_cells <- meta_data[order(meta_data$AUC, decreasing = TRUE), ][1:ceiling(nrow(meta_data) * 0.1), ]
  
  # Barplot for top 10% aging cells
  plot_top_cells_bar(top_cells, "top10_aging_cells_count.pdf")
  
  # Proportion barplot
  plot_aging_ratio(top_cells, meta_data, "aging_cell_proportion.pdf")
}
