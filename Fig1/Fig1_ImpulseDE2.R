rm(list = ls())

#----------------------------#
# Differential Gene Analysis with ImpulseDE2
#----------------------------#

# Load required packages
#devtools::install_github("YosefLab/ImpulseDE2")
library(ImpulseDE2)
library(stringr)
library(compiler)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)

# Set working directory
setwd("/home/data/t140311/Organ_aging")

# Load expression data
mydata <- readRDS("data/GSE132040_groupbysample.rds")

# Helper function to pad age values with leading zero
pad_age <- function(age_vec) {
  age_vec <- ifelse(nchar(age_vec) < 3, str_c("0", age_vec), age_vec)
  age_list <- strsplit(age_vec, "")
  return(as.numeric(sapply(age_list, function(x) paste0(x[1], x[2]))))
}


# Helper function to convert matrix to integer
force_matrix_to_integer <- function(m) {
  apply(m, c(1, 2), function(x) as.integer(x))
}

# Helper function to evaluate impulse model
eval_impulse <- cmpfun(function(vecImpulseParam, vecTimepoints) {
  values <- sapply(vecTimepoints, function(t) {
    (1 / vecImpulseParam[3]) * 
      (vecImpulseParam[2] + (vecImpulseParam[3] - vecImpulseParam[2]) /
         (1 + exp(-vecImpulseParam[1] * (t - vecImpulseParam[5])))) *
      (vecImpulseParam[4] + (vecImpulseParam[3] - vecImpulseParam[4]) /
         (1 + exp(vecImpulseParam[1] * (t - vecImpulseParam[6]))))
  })
  values[values < 1e-10] <- 1e-10
  return(values)
})

# Main loop for each organ
for (i in 1:16) {
  
  # Extract organ data
  data_organ <- as.data.frame(t(mydata[[i]]))
  organ <- data_organ[3, 1]
  
  # Create output directory
  result_dir <- file.path("/home/data/t140311/Organ_aging/result", organ)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(result_dir)
  
  # Prepare annotation data
  annotation <- data_organ[1:4,] %>% t() %>% as.data.frame()
  annotation$Condition <- "case"
  annotation$Time <- pad_age(as.vector(annotation$age))
  annotation$Batch <- sapply(strsplit(as.vector(t(data_organ[1, ])), "_"), `[`, 3)
  rownames(annotation) <- annotation$Sample
  colnames(annotation) <- c("Title","Sample", "organ","Age", "Time", "Condition")
  saveRDS(annotation, "Annotation_data.rds")
  
  # Process expression data
  exp_data <- data_organ[-c(2, 3, 4), ]
  colnames(exp_data) <- exp_data[1, ]
  exp_data <- exp_data[-1, ]
  rownames(exp_data) <- rownames(data_organ)[-c(1:4)]
  exp_data <- as.data.frame(lapply(exp_data, as.numeric))
  rownames(exp_data) <- rownames(data_organ)[-c(1:4)]
  exp_data <- force_matrix_to_integer(as.matrix(exp_data))
  saveRDS(as.data.frame(exp_data), "all_data_exp.rds")
  
  # Run ImpulseDE2
  objectImpulseDE2 <- runImpulseDE2(
    matCountData   = exp_data,
    dfAnnotation   = annotation,
    boolCaseCtrl   = FALSE,
    vecConfounders = NULL,
    scaNProc       = 1
  )
  saveRDS(objectImpulseDE2, "objectImpulseDE2.rds")
  
  # Extract DE results
  result <- objectImpulseDE2$dfImpulseDE2Results
  result_sig <- result[result$p < 0.01 & complete.cases(result), ]
  
  # Fit curves for DE genes
  gene_list <- rownames(result_sig)
  dat_line <- return_dfFit(
    vecGeneIDs       = gene_list[1],
    scaNTopIDs       = NULL,
    objectImpulseDE2 = objectImpulseDE2,
    boolCaseCtrl     = FALSE
  )[, 2, drop = FALSE]
  colnames(dat_line) <- gene_list[1]
  
  for (j in 2:length(gene_list)) {
    fit <- return_dfFit(
      vecGeneIDs       = gene_list[j],
      scaNTopIDs       = NULL,
      objectImpulseDE2 = objectImpulseDE2,
      boolCaseCtrl     = FALSE
    )
    rownames(fit) <- fit[, 1]
    fit <- fit[, 2, drop = FALSE]
    colnames(fit) <- gene_list[j]
    dat_line <- cbind(dat_line, fit)
  }
  
  # Filter genes with strictly monotonic trend
  Fit_data <- data.frame(test = rep(1, 100))
  for (j in 1:ncol(dat_line)) {
    diff_vec <- diff(as.numeric(dat_line[, j]))
    if (abs(sum(sign(diff_vec))) == 99) {
      Fit_data <- cbind(Fit_data, dat_line[, j, drop = FALSE])
    }
  }
  
  saveRDS(colnames(Fit_data)[-1], "diff_gene_filer.RDS")
  
  # Optional: Save DE result table
  write.table(result_sig, file = file.path(result_dir, "diff_genes.txt"), sep = "\t", row.names = FALSE)
  
  # Plot expression for selected genes (top 12)
  lsgplotsGenes <- plotGenes(
    vecGeneIDs       = colnames(Fit_data)[2:13],
    objectImpulseDE2 = objectImpulseDE2,
    boolCaseCtrl     = FALSE
  )
  print(lsgplotsGenes[[2]])
  
  p_combined <- cowplot::plot_grid(plotlist = lsgplotsGenes[1:12], nrow = 4, labels = LETTERS[1:12])
  print(p_combined)
  
  # Draw expression heatmap
  lsHeatmaps <- plotHeatmap(
    objectImpulseDE2       = objectImpulseDE2,
    strCondition           = "case",
    boolIdentifyTransients = FALSE,
    scaQThres              = 0.01
  )
  draw(lsHeatmaps$complexHeatmapRaw)
}
