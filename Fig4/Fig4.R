# Aging Score Analysis Pipeline with Trajectory-Based Gene Expression Evaluation

setwd("/data3/home/hu/Aging/Gemone medical")

# Load required libraries
library(Seurat)
library(monocle)
library(dplyr)
library(patchwork)
library(Matrix)
library(harmony)
library(GSEABase)
library(GSVA)
library(ggplot2)
library(cowplot)
library(ggridges)
library(ggpubr)
library(ggpmisc)
library(ggprism)
library(ggsci)

#------------------------------------------#
# Function: Create and preprocess Seurat object
#------------------------------------------#
prepare_seurat_object <- function(scRNA_path, gene_path, cell_number = 100000) {
  scRNA <- readRDS(scRNA_path)
  set.seed(123)
  if (ncol(scRNA) >= cell_number) {
    random_cells <- sample(colnames(scRNA), size = cell_number, replace = FALSE)
    scRNA_subset <- subset(scRNA, cells = random_cells)
  } else {
    stop("Insufficient number of cells.")
  }
  saveRDS(scRNA_subset, "Lung/scRNA_subset.RDS")
  
  gene <- read.csv(gene_path, header = TRUE)
  dup_gene <- gene[!duplicated(gene$gene_name), ]
  count_matrix <- scRNA_subset@assays$RNA$counts
  count_matrix <- count_matrix[!duplicated(gene$gene_name), ]
  rownames(count_matrix) <- dup_gene$gene_name
  
  metadata <- scRNA_subset@meta.data
  new_seurat <- CreateSeuratObject(counts = count_matrix, project = "New_Project")
  rownames(metadata) <- colnames(new_seurat)
  new_seurat <- AddMetaData(new_seurat, metadata = metadata)
  return(new_seurat)
}

#------------------------------------------#
# Function: Preprocess, integrate, and visualize data
#------------------------------------------#
process_and_plot <- function(seurat_object) {
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object)
  
  harmony_obj <- RunHarmony(seurat_object, "Age_group")
  harmony_obj <- RunUMAP(harmony_obj, reduction = "harmony", dims = 1:30)
  harmony_obj <- FindNeighbors(harmony_obj, reduction = "harmony", dims = 1:30)
  harmony_obj <- FindClusters(harmony_obj, resolution = 0.5)
  
  ncluster <- length(unique(harmony_obj$seurat_clusters))
  mycol <- colorRampPalette(brewer.pal(12, "Set3"))(ncluster)
  
  p1 <- DimPlot(harmony_obj, reduction = "umap", label = TRUE, pt.size = 0.5, cols = mycol, group.by = "seurat_clusters")
  p2 <- DimPlot(harmony_obj, reduction = "umap", group.by = "Main_cell_type", pt.size = 0.5)
  ggsave("umap_cluster.pdf", plot = p1 + p2, width = 10, height = 6)
  
  harmony_obj <- RunTSNE(harmony_obj, reduction = "harmony", dims = 1:30)
  harmony_obj <- FindNeighbors(harmony_obj, reduction = "harmony", dims = 1:30)
  harmony_obj <- FindClusters(harmony_obj, resolution = 0.5)
  
  p3 <- DimPlot(harmony_obj, reduction = "tsne", label = TRUE, pt.size = 0.5, cols = mycol, group.by = "seurat_clusters")
  p4 <- DimPlot(harmony_obj, reduction = "tsne", group.by = "Main_cell_type", pt.size = 0.5, cols = mycol)
  ggsave("tsne_cluster.pdf", plot = p3 + p4, width = 10, height = 6)
  
  return(harmony_obj)
}

#------------------------------------------#
# Function: Compute Aging Score
#------------------------------------------#
compute_aging_score <- function(scRNA_obj, up_gene, down_gene, organ_num_table) {
  exp_data <- t(as.data.frame(scRNA_obj@assays$RNA@layers$data))
  colnames(exp_data) <- rownames(scRNA_obj)
  
  up_gene_exp <- exp_data[, up_gene$gene]
  up_weights <- organ_num_table[organ_num_table$gene %in% up_gene$gene, ]$Organ_num * 0.1
  for (i in 1:ncol(up_gene_exp)) up_gene_exp[, i] <- up_gene_exp[, i] * up_weights[i]
  up_score <- rowMeans(up_gene_exp)
  
  down_gene_exp <- exp_data[, intersect(colnames(exp_data), down_gene$gene)]
  down_weights <- organ_num_table[organ_num_table$gene %in% down_gene$gene, ]$Organ_num * 0.1
  for (i in 1:ncol(down_gene_exp)) down_gene_exp[, i] <- down_gene_exp[, i] * down_weights[i]
  down_score <- rowMeans(down_gene_exp)
  
  age_score <- (up_score - down_score) / (sd(up_score - down_score) * sqrt(1/length(up_gene$gene) + 1/length(down_gene$gene)))
  return(age_score)
}

#------------------------------------------#
# Function: Visualize Aging Score
#------------------------------------------#
plot_aging_score <- function(meta_data, age_score) {
  meta_data$age_score <- age_score
  sample2 <- sapply(strsplit(as.character(meta_data$sample), "_"), function(x) x[7])
  meta_data$sample2 <- sample2
  meta_data$sample3 <- paste(meta_data$Age_group, sample2, sep = "_")
  meta_data$cell_type <- sapply(strsplit(as.character(meta_data$Sub_cell_type), "-"), function(x) x[1])
  
  cell_type_data <- meta_data %>%
    group_by(sample3, cell_type) %>%
    summarise(age_score = mean(age_score), .groups = "drop") %>%
    mutate(time_val = suppressWarnings(as.numeric(sub("_.*", "", sample3)))) %>%
    filter(!is.na(time_val)) %>%
    mutate(time = factor(time_val, levels = c(3, 6, 12, 16, 23)))
  
  
  p <- ggplot(cell_type_data, aes(x = time, y = age_score)) +
    stat_summary(fun = mean, geom = "bar", fill = "#858991", width = 0.7, alpha = 0.8) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, color = "black") +
    geom_smooth(aes(x = as.numeric(time), group = 1), method = "lm", formula = y ~ x, color = "#99BBE0", linetype = "dashed", se = FALSE, size = 0.8) +
    stat_cor(aes(x = as.numeric(time), y = age_score), method = "pearson", label.x.npc = "left", label.y.npc = 0.95, size = 4, color = "#99BBE0") +
    facet_wrap(~ cell_type, nrow = 4, scales = "free_y") +
    labs(x = "Age Group", y = "Aging Score") +
    theme_prism(base_line_size = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank())
  
  ggsave("Lung/age_group_barplot.pdf", plot = p, width = 16, height = 10)
}

new_obj <- prepare_seurat_object("GSE247719_PanSci_02_Lung_all.cell.RDS", "GSM7899265_PanSci_lung_df_gene.csv")
harmony_obj <- process_and_plot(new_obj)
age_score <- compute_aging_score(harmony_obj, up_gene, down_gene, group_Organ_num)
plot_aging_score(harmony_obj@meta.data, age_score)



######compare scImmuAging
library(scImmuAging)
# devtools::install_github("CiiM-Bioinformatics-group/scImmuAging")

model_set <- readRDS(system.file("data", "all_model.RDS", package = "scImmuAging"))
feature_set <- readRDS(system.file("data", "all_model_inputfeatures.RDS", package = "scImmuAging"))

feature_set1 <- list()
for(i in c("CD4T", "CD8T", "MONO", "NK", "B")) {
  temp_df <- coef(model_set[[i]])
  temp_feature <- rownames(temp_df)[which(temp_df[,1] != 0)]
  feature_set1[[i]] <- temp_feature
}    

all_features <- readRDS("Lung/aging_clock_comparison/all_features.RDS")
all_model <- readRDS("Lung/aging_clock_comparison/all_model.RDS")
all_model_input <- readRDS("Lung/aging_clock_comparison/all_model_inputfeatures.RDS")

preprocessing <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  meta_data <- seurat_obj@meta.data
  
  if(!("donor_id" %in% colnames(meta_data))) stop("donor_id missing in metadata!")
  if(!("age" %in% colnames(meta_data))) stop("age missing in metadata!")
  
  meta_data <- meta_data[, c("donor_id", "age")]
  input_mtx <- t(as.matrix(seurat_obj@assays$RNA$counts))
  combined_input <- as_tibble(cbind(meta_data, input_mtx))
  return(combined_input)
}

pseudocell <- function(input, size = 15, n = 100, replace = "dynamic") {
  pseudocells <- c()
  if (replace == "dynamic") {
    if (nrow(input) <= size) replace <- TRUE else replace <- FALSE
  }
  for (i in seq_len(n)) {
    batch <- input[sample(1:nrow(input), size = size, replace = replace), ]
    pseudocells <- rbind(pseudocells, colMeans(batch))
  }
  colnames(pseudocells) <- colnames(input)
  return(as_tibble(pseudocells))
}

count_matrix <- new_seurat@assays$RNA$counts

library(biomaRt)
library(dplyr)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

rownames(count_matrix) <- toupper(rownames(count_matrix))
gene_symbols <- rownames(count_matrix)

conversion_table <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = gene_symbols,
  mart = ensembl
)

conversion_table <- conversion_table[!duplicated(conversion_table$hgnc_symbol), ]
rownames(conversion_table) <- conversion_table$hgnc_symbol

common_genes <- intersect(conversion_table$hgnc_symbol, rownames(count_matrix))
count_matrix <- count_matrix[common_genes, ]
conversion_table <- conversion_table[common_genes, ]

rownames(count_matrix) <- conversion_table$ensembl_gene_id

metadata <- new_seurat@meta.data

sc_seurat <- CreateSeuratObject(
  counts = count_matrix,
  project = "New_Project"
)

rownames(metadata) <- colnames(sc_seurat)

new_seurat <- AddMetaData(
  object = sc_seurat,
  metadata = metadata
)

selected <- unique(c(all_features[[1]], all_features[[2]], all_features[[3]], all_features[[4]]))
sc_validation <- subset(sc_validation, features = intersect(selected, rownames(sc_validation)))

sc_validation$donor_id <- sapply(strsplit(as.character(sc_validation$sample), "_"), function(x) x[7])
sc_validation$donor_id <- paste(sc_validation$Age_group, sc_validation$donor_id, sep = "_")
sc_validation$age <- sc_validation$Age_group

internal_valid <- preprocessing(sc_validation) %>% group_by(donor_id, age) %>% nest()
internal_valid <- internal_valid %>% mutate(pseudocell_all = map(data, pseudocell))
internal_valid$data <- NULL
internal_valid <- unnest(internal_valid, pseudocell_all)

donor_id <- internal_valid[,1]
age <- internal_valid[,2]
final_mtx <- as.matrix(internal_valid[, -c(1,2)])

model <- all_model[["CD4T"]]
intersect_features <- intersect(colnames(final_mtx), rownames(model$glmnet.fit$beta))
nonintersect_features <- setdiff(colnames(final_mtx), rownames(model$glmnet.fit$beta))
model_len <- nrow(model$glmnet.fit$beta)
nonintersect_sampled <- sample(nonintersect_features, model_len - length(intersect_features))

final_mtx_temp <- final_mtx[, c(intersect_features, nonintersect_sampled)]
testPredictions <- predict(model, newx = final_mtx_temp, s = "lambda.min")

test_df <- data.frame(
  donor_id = donor_id,
  age = age,
  Prediction = testPredictions[,1]
)

test_df$age_num <- sapply(strsplit(as.character(test_df$age), "_"), function(x) x[1]) %>% as.numeric()

cor.test(test_df$age_num, test_df$Prediction)


###---------------------------------------------------##################


