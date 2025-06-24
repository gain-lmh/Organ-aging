#############################################
# Temporal Clustering using Mfuzz and Pathway Enrichment
#############################################

# Clear workspace and set working directory
rm(list = ls())
setwd("/home/data/t140311/Organ_aging")
organ_dirs <- dir()

# Example: process one organ directory
setwd(file.path("/home/data/t140311/Organ_aging", organ_dirs[11]))

# Load data and dependencies
library(ImpulseDE2)
library(limma)
library(edgeR)
library(ClusterGVis)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(gtools)
library(trend)
library(circlize)
library(ggsci)

# Step 1: Load data
exp_data <- readRDS("data_exp.rds")
anno_data <- readRDS("Annotation_data.rds")[, c("Sample", "Time")]
impulse_result <- readRDS("objectImpulseDE2.rds")

# Step 2: Select significantly DE genes
result <- impulse_result$dfImpulseDE2Results
result_filtered <- result[result$padj < 0.01 & complete.cases(result), ]
selected_exp <- exp_data[rownames(result_filtered), ]

# Step 3: Normalize and average expression by time
selected_exp <- as.data.frame(t(selected_exp))
selected_exp$Sample <- rownames(selected_exp)
merged_data <- merge(anno_data, selected_exp, by = "Sample")[, -1]
averaged_exp <- t(avereps(merged_data, ID = merged_data$Time))
colnames(averaged_exp) <- paste0("Month_", colnames(averaged_exp))
averaged_exp <- averaged_exp[-1, ]
averaged_exp <- averaged_exp[, paste0("Month_", c(1,3,6,9,12,15,18,21,24,27))]

# Step 4: Log2 CPM normalization
row_names <- rownames(averaged_exp)
averaged_exp <- log2(cpm(as.data.frame(lapply(averaged_exp, as.numeric))) + 1)
rownames(averaged_exp) <- row_names
saveRDS(averaged_exp, "avereps_df.rds")

# Step 5: Estimate optimal cluster number (manual check)
# getClusters(exp = averaged_exp)

# Step 6: Mfuzz clustering
c <- 12  ####Modify as needed
cm <- clusterData(exp = as.matrix(averaged_exp), cluster.method = "mfuzz", cluster.num = c)
membership_info <- data.frame(cm[[2]])
selected_genes <- membership_info[membership_info$membership > 0.5, ]
cm[[2]] <- selected_genes
exp_cluster <- data.frame(cm[[1]])
exp_cluster <- exp_cluster[exp_cluster$membership > 0.5, ]
cm[[1]] <- exp_cluster

# Save cluster gene list
gene_clusters <- selected_genes[, c("gene", "cluster")]
unique_genes <- gene_clusters[!duplicated(gene_clusters$gene), ]
write.table(unique_genes, paste0("cluster_", c, "_gene.txt"), sep = "\t", row.names = FALSE)

# Step 7: KEGG/GO Enrichment
gene_lists <- unique_genes %>% group_split(cluster)
names(gene_lists) <- paste0("C", 1:c)

convert_to_entrez <- function(gene_list) {
  s2e <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  return(s2e$ENTREZID)
}
gcSample <- lapply(gene_lists, function(df) convert_to_entrez(df$gene))

# KEGG enrichment
cm_KEGG <- compareCluster(gcSample, fun = "enrichKEGG", organism = "mmu", pvalueCutoff = 0.05)

# GO BP enrichment
cm_GO <- compareCluster(
  gcSample,
  fun = "enrichGO",
  OrgDb = "org.Mm.eg.db",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# Clean up KEGG descriptions
cm_KEGG <- as.data.frame(cm_KEGG)
cm_KEGG$Description <- sapply(strsplit(cm_KEGG$Description, " - Mus"), `[`, 1)

# Save top 5 terms per cluster
termanno_kegg <- cm_KEGG %>% group_by(Cluster) %>% top_n(-5, wt = pvalue) %>% select(Cluster, Description)
termanno_go <- as.data.frame(cm_GO) %>% group_by(Cluster) %>% top_n(-5, wt = pvalue) %>% select(Cluster, Description)
saveRDS(termanno_kegg, "termanno_kegg.rds")
saveRDS(termanno_go, "termanno_go.rds")

# Step 8: Cluster visualization
saveRDS(cm, paste0("cluster_", c, "_cm.rds"))

pdf(paste0("cluster_", c, "_gene.pdf"), height = 6, width = 8)
visCluster(object = cm, plot.type = "line", ncol = 3)
dev.off()

colormap <- colorRamp2(breaks = c(-2.5, 0, 2.5), colors = c("#328785", "white", "#C78177"))
pdf(paste0("all_cluster_test_", c, "_gene.pdf"), height = 10, width = 12)
visCluster_plot(
  object = cm,
  plot.type = "both",
  column_names_rot = 45,
  annoTerm.data = termanno_kegg,
  line.side = "left",
  show_row_dend = FALSE,
  heatmap.col = colormap,
  add.box = TRUE,
  boxcol = pal_npg()(8)
)
dev.off()

# Step 9: Mann-Kendall test for monotonic trends
mk_pvals <- mk_stats <- c()
for (i in 1:length(gene_lists)) {
  clust_data <- exp_cluster[exp_cluster$cluster == i, ]
  max_row <- which.max(clust_data$membership)
  expression <- t(clust_data[max_row, 2:11])
  mk_result <- mk.test(expression[, 1], continuity = TRUE)
  mk_pvals <- c(mk_pvals, mk_result$p.value)
  mk_stats <- c(mk_stats, mk_result$statistic)
}

significant_clusters <- which(mk_pvals <= 0.05)
trend_genes <- exp_cluster[exp_cluster$cluster %in% significant_clusters, ]
trend_genes$membership <- trend_genes$membership / ave(trend_genes$membership, trend_genes$cluster, FUN = max)
filtered_genes <- trend_genes[trend_genes$membership > 0.8, ]
cm$wide.res <- filtered_genes

pdf(paste0("HTtermCmlsrt_cluster_", c, "_gene.pdf"), height = 8, width = 12)
visCluster(
  object = cm,
  plot.type = "both",
  column_names_rot = 45,
  annoTerm.data = termanno_kegg,
  line.side = "left",
  show_row_dend = FALSE
)
dev.off()

filtered_genes$mk.statistic <- rep(mk_stats[significant_clusters], times = table(filtered_genes$cluster))
write.table(filtered_genes, paste0("cluster_", c, "_file_gene_exp.txt"))
saveRDS(filtered_genes, paste0("cluster_", c, "_mfuzz_file_gene_exp.rds"))
saveRDS(cm_GO, "cm_GO.rds")
