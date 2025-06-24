rm(list = ls())

# Load required packages
suppressMessages({
  library(devtools)
  library(BayesPrism)
  library(Seurat)
  library(ggplot2)
  library(ComplexHeatmap)
  library(Hmisc)
})

# Set file paths
base_path <- "/home/data/t160407/Organ_aging/Mouse Aging/Spleen"
bulk_file <- file.path(base_path, "Row data/Bulk/spleen_exp_GEO.txt")
scRNA_file <- file.path(base_path, "Result/spleen_after_anno_scRNA_harmony.RDS")

# Load bulk data
bulk_data <- read.table(bulk_file, header = TRUE, sep = "\t")
bulk_data <- bulk_data[!is.na(bulk_data$GeneSymbol), ]
rownames(bulk_data) <- bulk_data$GeneSymbol
bulk_data <- t(bulk_data[, -1])
bulk_data <- as.data.frame(bulk_data)

# Load single-cell data
sc_data <- readRDS(scRNA_file)
Idents(sc_data) <- sc_data$category
cell_labels <- sc_data$category
sc_matrix <- t(as.matrix(sc_data@assays$RNA@counts))
sc_matrix <- as.data.frame(sc_matrix)

# Visualize cell type correlation
plot.cor.phi(
  input = sc_matrix,
  input.labels = cell_labels,
  title = "Cell Type Correlation"
)

# Visualize outlier genes in scRNA and bulk
sc_outlier <- plot.scRNA.outlier(
  input = sc_matrix,
  cell.type.labels = cell_labels,
  species = "mm",
  return.raw = TRUE
)

bulk_outlier <- plot.bulk.outlier(
  bulk.input = bulk_data,
  sc.input = sc_matrix,
  cell.type.labels = cell_labels,
  species = "mm",
  return.raw = TRUE
)

# Filter outlier genes in single-cell data
sc_matrix_filtered <- cleanup.genes(
  input = sc_matrix,
  input.type = "count.matrix",
  species = "mm",
  gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
  exp.cells = 5
)

# Check consistency between bulk and scRNA expression
plot.bulk.vs.sc(
  sc.input = sc_matrix_filtered,
  bulk.input = bulk_data
)

# Extract protein-coding genes only
sc_matrix_pc <- select.gene.type(
  sc_matrix_filtered,
  gene.type = "protein_coding"
)

# Build BayesPrism model
prism_model <- new.prism(
  reference = sc_matrix_pc,
  mixture = bulk_data,
  input.type = "count.matrix",
  cell.type.labels = cell_labels,
  cell.state.labels = NULL,
  key = NULL,
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

# Run BayesPrism inference
prism_result <- run.prism(prism = prism_model, n.cores = 50)

# Extract estimated cell-type proportions
theta <- get.fraction(
  bp = prism_result,
  which.theta = "final",
  state.or.type = "type"
)

theta <- as.data.frame(theta)

# Define sample time groups (example: 6 time points × 3 replicates)
time_group <- rep(c("13", "26", "52", "78", "104", "130"), each = 3)

# Correlation analysis: time vs. cell-type proportion
cor_plots <- list()
p_values <- c()

for (i in seq_len(ncol(theta))) {
  temp_df <- data.frame(
    proportion = theta[, i],
    group = factor(time_group, levels = c("13", "26", "52", "78", "104", "130"))
  )
  temp_df$group_numeric <- as.numeric(temp_df$group)
  
  cor_result <- rcorr(temp_df$group_numeric, temp_df$proportion)
  p <- cor_result$P[1, 2]
  r <- round(cor_result$r[1, 2], 2)
  p_text <- ifelse(p < 0.01, "p < 0.01", paste0("p = ", signif(p, 2)))
  
  plot <- ggplot(temp_df, aes(x = group_numeric, y = proportion)) +
    geom_point(color = "#988d7b") +
    geom_smooth(method = "lm", formula = y ~ x, fill = "#b2e7fa", color = "#00aeef", alpha = 0.8) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold.italic"),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = paste0(colnames(theta)[i], "  ρ = ", r, ", ", p_text))
  
  cor_plots[[i]] <- plot
  p_values <- c(p_values, p)
}

# Display significant cell types
significant_types <- colnames(theta)[p_values < 0.05]
print("Significant cell types:")
print(significant_types)

# Plot heatmap of inferred proportions
Heatmap(theta, name = "Proportion", row_km = 3, column_km = 3)
