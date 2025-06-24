############################################
# Consensus Clustering
############################################

rm(list = ls())
library(ConsensusClusterPlus)
library(networkD3)

#----------------------------#
# Function: Load and Prepare Expression Data
#----------------------------#
load_expression_data <- function(expr_path, gene_list) {
  expr_data <- readRDS(expr_path)
  filtered_expr <- expr_data[gene_list, ]
  log2_expr <- log2(filtered_expr + 1)
  return(as.matrix(log2_expr))
}

#----------------------------#
# Function: Perform Consensus Clustering
#----------------------------#
run_consensus_clustering <- function(expr_matrix, k = 4, seed = 1234) {
  expr_centered <- sweep(expr_matrix, 1, apply(expr_matrix, 1, median))
  
  cluster_result <- ConsensusClusterPlus(
    d = expr_centered,
    maxK = k,
    pItem = 0.8,
    pFeature = 1,
    clusterAlg = "hc",
    distance = "pearson",
    seed = seed,
    innerLinkage = "complete",
    finalLinkage = "complete",
    corUse = "pairwise.complete.obs",
    plot = "pdf",
    title = "consensus_clustering"
  )
  
  return(cluster_result)
}

#----------------------------#
# Function: Create Cluster Assignment Metadata
#----------------------------#
create_cluster_metadata <- function(anno_data, cluster_assignment) {
  anno_data$Time <- paste0(anno_data$Time, "_month")
  metadata <- data.frame(
    sample = anno_data$Sample,
    Time = anno_data$Time,
    group = cluster_assignment,
    value = 1
  )
  return(metadata)
}

#----------------------------#
# Function: Prepare Sankey Node and Link Data
#----------------------------#
prepare_sankey_data <- function(cluster_data) {
  cluster_data$group <- as.factor(cluster_data$group)
  cluster_data$Time <- as.factor(cluster_data$Time)
  
  factors_all <- sort(unique(c(levels(cluster_data$Time), levels(cluster_data$group))))
  node_ids <- 0:(length(factors_all) - 1)
  
  levels(cluster_data$Time) <- node_ids[factors_all %in% levels(cluster_data$Time)]
  levels(cluster_data$group) <- node_ids[factors_all %in% levels(cluster_data$group)]
  
  cluster_data$Time <- as.numeric(as.character(cluster_data$Time))
  cluster_data$group <- as.numeric(as.character(cluster_data$group))
  
  node_df <- data.frame(name = factors_all)
  return(list(links = cluster_data, nodes = node_df))
}

#----------------------------#
# Function: Plot Sankey Diagram
#----------------------------#
plot_sankey <- function(links, nodes) {
  sankey_color <- 'd3.scaleOrdinal()
    .domain(["1", "2", "3", "4", "6", "7", "8", "9", "10", "11", "0", "5"])
    .range(["#54686F", "#E57B7F", "#CB997E", "#d0b9d5", "#9E3150", "#B8B7A3",
            "#60b4e2", "#FC8C5A", "#D44C3C", "#B1EAEF", "#4EAB90", "#8EB69C"])'
  
  p <- sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "Time",
    Target = "group",
    Value = "value",
    NodeID = "name",
    fontSize = 12,
    nodeWidth = 30,
    colourScale = sankey_color
  )
  return(p)
}

#----------------------------#
# MAIN PIPELINE EXECUTION
#----------------------------#

# Set working directory
setwd("/home/data/t140311/Organ_aging")
file_list <- list.files(full.names = TRUE)
filtered_gene_file <- readRDS(file_list[2])

# Load organ-specific expression and annotation data
setwd("/home/data/t140311/Organ_aging/Brain/")
expression_matrix <- load_expression_data("data_exp.rds", filtered_gene_file$gene)
annotation_data <- readRDS("Annotation_data.rds")

# Run clustering
cluster_result <- run_consensus_clustering(expression_matrix, k = 4)
cluster_assignment <- cluster_result[[3]]$consensusClass

# Generate metadata
cluster_metadata <- create_cluster_metadata(annotation_data, cluster_assignment)
saveRDS(cluster_metadata, "consensus_cluster.rds")

# Prepare data and plot
sankey_components <- prepare_sankey_data(cluster_metadata)
sankey_plot <- plot_sankey(sankey_components$links, sankey_components$nodes)

