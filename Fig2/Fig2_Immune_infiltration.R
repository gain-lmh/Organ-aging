rm(list = ls())

library(edgeR)
library(CIBERSORT)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(ggthemes)
library(jjPlot)
library(jjAnno)
library(readxl)
library(ggsci)

set_base_dirs <- function() {
  base_dir <- "D:/result/"
  mouse_ref <- "D:/data/share/mice.txt"
  list(base_dir = base_dir, mouse_ref = mouse_ref)
}

#------------------------------#
# Step 1: Preprocess Expression
#------------------------------#
normalize_expression <- function(expr_matrix) {
  rownames(expr_matrix) <- rownames(expr_matrix)
  expr_matrix <- expr_matrix[rowSums(expr_matrix >= 1) >= 50, colSums(expr_matrix >= 1) >= 10000]
  expr_matrix <- log2(cpm(as.matrix(expr_matrix)) + 1)
  return(expr_matrix)
}

#------------------------------#
# Step 2: Run CIBERSORT
#------------------------------#
run_cibersort_analysis <- function(expr_matrix, signature_matrix) {
  cibersort(sig_matrix = as.matrix(signature_matrix),
            mixture_file = as.matrix(expr_matrix),
            perm = 1000,
            QN = TRUE)
}

#------------------------------#
# Step 3: Normalize Cell Fractions (row-wise)
#------------------------------#
normalize_cibersort_output <- function(ciber_result) {
  result <- ciber_result[, 1:(ncol(ciber_result)-3)]  # Exclude P-value, Correlation, RMSE
  result <- result[, colSums(result == 0) < 30]
  for (i in 1:ncol(result)) {
    result[, i] <- result[, i] / sum(result[, i])
  }
  return(result)
}

#------------------------------#
# Step 4: Extract and Save
#------------------------------#
process_organ <- function(folder, signature_matrix) {
  setwd(folder)
  cluster_data <- readRDS("consensus_cluster.rds")
  exp_data <- readRDS("data_exp.rds")
  norm_expr <- normalize_expression(exp_data)
  ciber_res <- run_cibersort_analysis(norm_expr, signature_matrix)
  cibersort_norm <- normalize_cibersort_output(ciber_res)
  saveRDS(cibersort_norm, "mouse_immune_infiltration/immune_profile.rds")
  return(cibersort_norm)
}

#------------------------------#
# Step 5: Collect All Organs
#------------------------------#
run_all_organs <- function(base_dir, signature_matrix) {
  folders <- list.dirs(base_dir, recursive = FALSE)
  immune_list <- list()
  
  for (folder in folders) {
    setwd(folder)
    dir.create("mouse_immune_infiltration", showWarnings = FALSE)
    message("Processing: ", folder)
    immune_list[[basename(folder)]] <- process_organ(folder, signature_matrix)
  }
  return(immune_list)
}

#------------------------------#
# Step 6: Statistical Testing
#------------------------------#
compare_groups_wilcox <- function(cluster_df, immune_df, organ_name) {
  merged_df <- merge(cluster_df[, c("sample", "group")], immune_df, by = "sample")
  merged_df <- merged_df[, colSums(merged_df == 0) < 30]
  
  pvals <- sapply(3:ncol(merged_df), function(j) {
    wilcox.test(merged_df[merged_df$group == 1, j], 
                merged_df[merged_df$group == 2, j])$p.value
  })
  
  group_avg <- merged_df %>%
    select(-sample) %>%
    group_by(group) %>%
    summarise(across(everything(), mean)) %>%
    mutate(Group = ifelse(group == 1, paste0(organ_name, "_young"), paste0(organ_name, "_old"))) %>%
    select(-group)
  
  mat <- as.data.frame(t(group_avg[, -ncol(group_avg)]))
  mat$immune_type <- rownames(mat)
  mat$p.value <- pvals
  mat$organ <- organ_name
  return(mat)
}

#------------------------------#
# Step 7: Merge Results
#------------------------------#
merge_result_lists <- function(result_list) {
  Reduce(function(x, y) merge(x, y, all = TRUE), result_list)
}

#------------------------------#
# Step 8: Visualize Significant Differences
#------------------------------#
plot_significance_bar <- function(data, axis_label, file_name) {
  ggplot(data, aes(x = reorder(type, x), y = x)) +
    geom_bar(stat = "identity", fill = "#479E9B", width = 0.7) +
    geom_text(aes(label = x), vjust = -0.5, size = 3.5) +
    labs(y = "Significant Immune Sets", x = axis_label) +
    coord_flip() +
    theme_classic() +
    theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1)) +
    ggsave(file_name, width = 6, height = 4)
}

#------------------------------#
# MAIN EXECUTION
#------------------------------#

dirs <- set_base_dirs()
mouse_signature <- read.table(dirs$mouse_ref, sep = "\t", header = TRUE, row.names = 1, fill = TRUE, fileEncoding = "UTF-16LE")
immune_list <- run_all_organs(dirs$base_dir, mouse_signature)

# Load processed CIBERSORT results for significance testing
organ_dirs <- list.dirs(dirs$base_dir, recursive = FALSE)
pval_results <- list()

for (org in organ_dirs) {
  organ_name <- basename(org)
  setwd(org)
  cluster <- readRDS("consensus_cluster.rds")
  immune <- readRDS("mouse_immune_infiltration/immune_profile.rds")
  immune$sample <- rownames(immune)
  res <- compare_groups_wilcox(cluster, immune, organ_name)
  pval_results[[organ_name]] <- res
}

# Merge all p-value results
all_pval <- do.call(rbind, pval_results)

# Count significant immune cell types by organ and by immune type
all_pval$significant <- ifelse(all_pval$p.value < 0.05, 1, 0)

sig_by_organ <- all_pval %>%
  group_by(organ) %>%
  summarise(sig_count = sum(significant)) %>%
  rename(type = organ, x = sig_count)

sig_by_immune <- all_pval %>%
  group_by(immune_type) %>%
  summarise(sig_count = sum(significant)) %>%
  rename(type = immune_type, x = sig_count)

# Plot summary bar charts
plot_significance_bar(sig_by_organ, "Organ", "organ_significance.pdf")
plot_significance_bar(sig_by_immune, "Immune Cell Type", "immune_type_significance.pdf")
