rm(list = ls())
### Step1 Data Preprocessing ####
{
  ### Load drug-target mapping data from DrugBank database
  Drug_bank <- read.csv("total data/drug links.csv")  #### Mapping between drug IDs and drug names
  
  ### Load annotation info from CMap database
  sig_anno <- read.table("total data/siginfo_beta.txt", header = TRUE,
                         sep = "\t") ##### Annotation info in CMap database
  
  ### Load cell line annotation from CMap database
  celltype_anno <- read.csv("total data/celltype.csv") ##### Cell line annotations
  
  ### Select drugs that have gene expression perturbation data
  Drugbank_Camp_interDrug <- intersect(capitalize(sig_anno$cmap_name), Drug_bank$Name)
  
  #### Extract intersected drug gene expression data — too large, load on demand
  gctx_demo = parse_gctx("level5_beta_trt_cp_n720216x12328.gctx")
  gctx_demo@mat[1:4,1:4]
  
  ### Drug intersection annotation info
  temp_sig_anno <- sig_anno[capitalize(sig_anno$cmap_name) %in% Drugbank_Camp_interDrug, ]
  temp_sig_anno_24h <- temp_sig_anno[temp_sig_anno$pert_time %in% "24", ]
  
  ### Select normal cell lines
  normal_celltype <- celltype_anno[celltype_anno$cell_type == "normal", ]
  
  #### Filter drugs with 24-hour perturbation data in normal cell lines
  temp_sig_anno_24h_normol <- temp_sig_anno_24h[temp_sig_anno_24h$cell_iname %in% normal_celltype$cell_iname, ]
  
  #### Control group gene expression data
  Ctrl_exp <- parse_gctx("total data/level5_beta_ctl_n58022x12328.gctx")
  Ctrl_exp@mat[1:4,1:4]
  Ctrl_exp_name <- Ctrl_exp@cid
  
  #### Control group annotation info for intersected drugs, restricted to normal cell lines
  Ctrl_sig_anno <- sig_anno[capitalize(sig_anno$sig_id) %in% Ctrl_exp_name, ]
  Ctrl_sig_anno_normal <- Ctrl_sig_anno[Ctrl_sig_anno$cell_iname %in% normal_celltype$cell_iname, ]
  
  #### Load drug-target data from DrugBank
  Drug_bank_target <- read.csv("total data/drug_target_final.csv")  #### Drug-target dataset
  Drug_bank_total <- merge(Drug_bank_target, Drug_bank, by.x = "drug_ID_1", by.y = "DrugBank.ID")
  
  #### Re-filter drugs intersecting DrugBank and CMap 24h data in normal cell lines
  Drugbank_Camp_interDrug_24h <- intersect(capitalize(temp_sig_anno_24h_normol$cmap_name), Drug_bank$Name)
  
  Drug_Camp_bank_target <- Drug_bank_total[Drug_bank_total$Name %in% Drugbank_Camp_interDrug_24h, ]
  Drug_Camp_bank_target <- Drug_Camp_bank_target[, c(1,3,4,5,6)]
  
  #### Save processed data for later steps
  saveRDS(Drug_Camp_bank_target, "Step1 data/Drug_Camp_bank_target.RDS") #### Drug-target intersection data
  saveRDS(temp_sig_anno_24h_normol, "Step1 data/temp_sig_anno_24h_normol.RDS") #### Drug perturbation annotation in CMap (normal cell lines, 24h)
  saveRDS(Ctrl_sig_anno_normal, "Step1 data/Ctrl_sig_anno_normal.RDS") #### Control group annotation in CMap (normal cell lines)
}


#### Step2 Random Walk Scoring of Drugs Based on Aging Trend Genes ####
{
  setwd("/home/data/t160407/Organ_aging/Drug Filter/Aging Gene/")
  rm(list = ls())
  dir <- dir()
  library(RandomWalkRestartMH)
  library(igraph)
  mol_result_list_score <- list()
  
  for (i in 1:length(dir)) {
    
    #### Load hypergeometric test filtered drugs
    Drug_Camp_bank_target <- readRDS("..//Step1 data/Drug_Camp_bank_target.RDS")
    
    #### Load and process PPI network
    PPI_file_name <- paste0(dir[[i]], "//", dir[[i]], ".tsv")
    
    library(data.table)
    String_data <- fread(file = PPI_file_name, sep = '\t', header = TRUE, check.names = FALSE)
    String_data <- String_data[, c(1, 2)]
    
    point_data <- unique(c(String_data$`#node1`, String_data$node2))
    string_net <- graph_from_data_frame(d = String_data, vertices = point_data, directed = TRUE) 
    
    #### Calculate page rank values of aging trend genes in PPI network
    page_rank <- page.rank(string_net, algo = c("prpack"), 
                           vids = V(string_net), directed = TRUE, damping = 0.85)
    page_rank$vector <- page_rank$vector * 100000
    
    #### Build drug-target molecular network
    molecule_target <- Drug_Camp_bank_target[, c(4, 3)]
    colnames(molecule_target) <- c("node1", "node2")
    colnames(String_data) <- c("node1", "node2")
    
    molecule_target_string <- rbind(molecule_target, String_data[, c(1, 2)])  ## Integrated small molecule-target network
    
    #### Construct igraph object for multiplex network
    point_data <- unique(c(molecule_target_string$node2, molecule_target_string$node1))
    net <- graph_from_data_frame(d = molecule_target_string, vertices = point_data, directed = TRUE)
    
    #### Create Multiplex object and compute normalized adjacency matrix
    PPI_MultiplexObject <- create.multiplex(list(PPI = net))
    
    AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
    AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)
    
    #### Custom function to generate seed scores for multiplex random walk
    get.seed.scoresMultiplex <- function(Seeds, Number_Layers, tau) {
      Nr_Seeds <- length(Seeds)
      Seeds_Seeds_Scores <- rep(tau / Nr_Seeds, Nr_Seeds)
      Seed_Seeds_Layer_Labeled <- paste0(rep(Seeds, Number_Layers), sep = "_", rep(seq(Number_Layers), length.out = Nr_Seeds * Number_Layers, each = Nr_Seeds))
      Seeds_Score <- data.frame(Seeds_ID = Seed_Seeds_Layer_Labeled, Score = Seeds_Seeds_Scores, stringsAsFactors = FALSE)
      return(Seeds_Score)
    }
    
    #### Function to compute geometric mean of scores over layers
    geometric.mean <- function(Scores, L, N) {
      FinalScore <- numeric(length = N)
      for (i in seq_len(N)) {
        FinalScore[i] <- prod(Scores[seq(from = i, to = N * L, by = N)])^(1 / L)
      }
      return(FinalScore)
    }
    
    #### Parameters for random walk
    x <- AdjMatrixNorm_PPI
    MultiplexObject <- PPI_MultiplexObject
    L <- MultiplexObject$Number_of_Layers
    N <- MultiplexObject$Number_of_Nodes
    r <- 0.7
    tau <- rep(1, L) / L
    MeanType <- "Geometric"
    DispResults <- "TopScores"
    
    Threshold <- 1e-10
    NetworkSize <- ncol(x)
    residue <- 1
    iter <- 1
    
    #### Initialize seed scores and restart vector for random walk
    Seeds <- names(page_rank$vector)
    Seeds_Score <- data.frame(Seeds_ID = names(page_rank$vector), Score = page_rank$vector)
    
    prox_vector <- matrix(0, nrow = NetworkSize, ncol = 1)
    colname <- colnames(x) %>% strsplit("_") %>% sapply(function(x) x[1])
    prox_vector[which(colname %in% Seeds_Score[,1])] <- Seeds_Score[,2]
    
    prox_vector <- prox_vector / sum(prox_vector)
    restart_vector <- prox_vector
    
    #### Perform random walk with restart until convergence
    while(residue >= Threshold){
      old_prox_vector <- prox_vector
      prox_vector <- (1 - r) * (x %*% prox_vector) + r * restart_vector
      residue <- sqrt(sum((prox_vector - old_prox_vector)^2))
      iter <- iter + 1
    }
    
    #### Aggregate scores across layers by geometric mean
    rank_global <- data.frame(NodeNames = character(length = N), Score = numeric(length = N))
    rank_global$NodeNames <- gsub("_1", "", row.names(prox_vector)[seq_len(N)])
    
    if (MeanType == "Geometric"){
      rank_global$Score <- geometric.mean(as.vector(prox_vector[,1]), L, N)
    } else if (MeanType == "Arithmetic") {
      rank_global$Score <- regular.mean(as.vector(prox_vector[,1]), L, N)
    } else {
      rank_global$Score <- sumValues(as.vector(prox_vector[,1]), L, N)
    }
    
    #### Sort and filter results
    if (DispResults == "TopScores"){
      Global_results <- rank_global[with(rank_global, order(-Score, NodeNames)), ]
      Global_results <- Global_results[which(!Global_results$NodeNames %in% Seeds),]
    } else {
      Global_results <- rank_global    
    }
    
    rownames(Global_results) <- c()
    
    RWRM_ranking <- list(RWRM_Results = Global_results, Seed_Nodes = Seeds)
    class(RWRM_ranking) <- "RWRM_Results"
    
    mol_result <- RWRM_ranking$RWRM_Results
    mol_result <- mol_result[mol_result$NodeNames %in% Drug_Camp_bank_target$Name, ]
    mol_result <- mol_result[mol_result$Score > max(mol_result$Score) / 500, ]
    
    mol_result_list_score[[i]] <- mol_result
  }
  
  names(mol_result_list_score) <- dir
  saveRDS(mol_result_list_score, "../Step2 data/mol_result_list_score.RDS")
}

#### Step3 GSEA Enrichment Analysis to Obtain Drug Scores ####
rm(list = ls())
#### Load necessary data
setwd("/home/data/t160407/Organ_aging/Drug Filter")
Trend_gene_list <- readRDS("Step1 data/Trend_gene_list.RDS") ### Aging trend gene list

Trend_gene_list <- Trend_gene_list[-11]
sig_anno_24h_normol <- readRDS("Step1 data/temp_sig_anno_24h_normol.RDS") ### Drug perturbation annotation info
sig_anno_24h_normol$cmap_name <- capitalize(sig_anno_24h_normol$cmap_name)

Ctrl_sig_anno_normal <- readRDS("Step1 data/Ctrl_sig_anno_normal.RDS") ### Control group drug annotation info
mol_result_list_score <- readRDS("Step2 data/mol_result_list_score.RDS") #### Drug scores from Step2
mol_result_list_score <- mol_result_list_score[-4]

{
  library("clusterProfiler")
  
  gene <- rownames(merge_exp)
  gene.df <- bitr(gene, fromType = "ENTREZID", 
                  toType = c("SYMBOL"),
                  OrgDb = org.Hs.eg.db)
  head(gene.df)
  write.table(gene.df, "total data/gene.df.txt", sep = "\t")
  
  Organ_drug_Gsea <- list()
  
  for (i in 1:length(mol_result_list_score)) {
    #### Extract drug names for this organ
    Organ_drug <- mol_result_list_score[[i]]$NodeNames
    temp_anno <- sig_anno_24h_normol[sig_anno_24h_normol$cmap_name %in% Organ_drug, ]
    Gsea_NES_list <- list()
    
    for (j in 1:length(Organ_drug)) {
      ### Drug perturbation samples
      One_drug_anno <- temp_anno[temp_anno$cmap_name %in% Organ_drug[j], ]
      max_rep <- names(table(One_drug_anno$nearest_dose))[table(One_drug_anno$nearest_dose) == max(table(One_drug_anno$nearest_dose))]
      One_drug_anno <- One_drug_anno[One_drug_anno$nearest_dose == max_rep[length(max_rep)], ]
      
      ### Control group samples from same cell lines and 24h time point
      Ctrl_anno <- Ctrl_sig_anno_normal[Ctrl_sig_anno_normal$cell_iname %in% unique(One_drug_anno$cell_iname), ]
      Ctrl_anno <- Ctrl_anno[Ctrl_anno$pert_itime == "24 h", ]
      
      #### Extract expression data of drug treatment group
      gctx_demo = parse_gctx("total data/level5_beta_trt_cp_n720216x12328.gctx", cid = One_drug_anno$sig_id)
      gctx_demo_exp <- gctx_demo@mat
      
      #### Extract expression data of control group
      Ctrl_demo = parse_gctx("total data/level5_beta_ctl_n58022x12328.gctx", cid = Ctrl_anno$sig_id)
      Ctrl_demo_exp <- Ctrl_demo@mat
      
      ### Merge treated and control expression data
      merge_exp <- cbind(gctx_demo_exp, Ctrl_demo_exp)
      merge_exp <- as.data.frame(merge_exp)
      
      #### Map gene IDs from ENTREZID to SYMBOL
      gene_id <- read.csv("total data/gene.df.csv", header = TRUE, row.names = 1)
      
      merge_exp$ENTREZID <- rownames(merge_exp)
      gene_merge_exp <- merge(gene_id, merge_exp, by = "ENTREZID")
      gene_merge_exp <- gene_merge_exp[!duplicated(gene_merge_exp$SYMBOL), ]
      rowname <- gene_merge_exp$SYMBOL
      gene_merge_exp <- gene_merge_exp[, c(-1, -2)]
      gene_merge_exp <- as.data.frame(lapply(gene_merge_exp, as.numeric))
      rownames(gene_merge_exp) <- rowname
      
      #### Define sample groups for differential expression
      Group <- c(rep("Drug", ncol(gctx_demo_exp)), rep("Ctrl", ncol(Ctrl_demo_exp))) %>% 
        factor(., levels = c("Drug", "Ctrl"), ordered = FALSE)
      
      #### Perform differential expression analysis using limma
      library(limma)
      library(dplyr)
      
      Group <- model.matrix(~factor(Group) + 0)  # Create model matrix
      colnames(Group) <- c("Drug", "Ctrl")
      df.fit <- lmFit(gene_merge_exp, Group)  # Fit linear model
      
      df.matrix <- makeContrasts(Drug - Ctrl, levels = Group)
      fit <- contrasts.fit(df.fit, df.matrix)
      fit <- eBayes(fit)
      tempOutput <- topTable(fit, n = Inf, adjust = "fdr")
      
      #### Extract aging trend genes up and downregulated
      Trend_gene <- Trend_gene_list[[i]]
      Up_Trend_gene <- Trend_gene[Trend_gene$mk.statistic > 0, ]$HGNC.symbol
      Down_Trend_gene <- Trend_gene[Trend_gene$mk.statistic < 0, ]$HGNC.symbol
      
      ### GSEA enrichment analysis
      gene_list <- tempOutput$logFC
      names(gene_list) <- rownames(tempOutput)
      gene_list <- sort(gene_list, decreasing = TRUE)
      
      TERM2GENE = data.frame(
        Term = c(rep("Rising", length(Up_Trend_gene)), rep("Declining", length(Down_Trend_gene))), 
        Gene = c(Up_Trend_gene, Down_Trend_gene)
      )
      
      gsea_result <- GSEA(geneList = gene_list, TERM2GENE = TERM2GENE,
                          pvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, verbose = FALSE)
      
      gsea_result <- gsea_result@result[, c(1,5,6,7)]
      gsea_result$drug <- rep(Organ_drug[j], nrow(gsea_result))
      Gsea_NES_list[[j]] <- gsea_result
    }
    
    names(Gsea_NES_list) <- Organ_drug
    combined_Gsea_NES <- do.call(rbind, Gsea_NES_list)
    Organ_drug_Gsea[[i]] <- combined_Gsea_NES
  }
}
names(Organ_drug_Gsea) <- names(mol_result_list_score)
saveRDS(Organ_drug_Gsea, "Step2 data/Organ_drug_Gsea.RDS")


##### Filter drugs that meet anti-aging criteria based on GSEA results #####
setwd("/home/data/t160407/Organ_aging/Drug Filter/Aging Gene")
dir <- dir()
dir <- dir[-4]
after_filter_GSEA_list <- list()

for(i in 1:length(Organ_drug_Gsea)){
  One_drug_Gsea <- Organ_drug_Gsea[[i]]
  
  ### Drugs with negative enrichment for Rising genes and positive for Declining genes are considered anti-aging
  negative_Rising <- One_drug_Gsea[One_drug_Gsea$ID == "Rising" & One_drug_Gsea$NES < 0, ]
  positive_Declining <- One_drug_Gsea[One_drug_Gsea$ID == "Declining" & One_drug_Gsea$NES > 0, ]
  
  anti_aging_Drug <- union(negative_Rising$drug, positive_Declining$drug)
  
  ### Add anti-aging or pro-aging label
  One_drug_Gsea$aging_anti <- ifelse(One_drug_Gsea$drug %in% anti_aging_Drug, "anti_aging", "posi_aging")
  One_drug_Gsea$abs_NES <- abs(One_drug_Gsea$NES)
  
  ### Sort by absolute NES descending
  One_drug_Gsea <- One_drug_Gsea[order(One_drug_Gsea$abs_NES, decreasing = TRUE), ]
  One_drug_Gsea <- One_drug_Gsea[!duplicated(One_drug_Gsea$drug), ]
  
  ### Add organ and rank info
  One_drug_Gsea$Organ <- rep(dir[i], nrow(One_drug_Gsea))
  One_drug_Gsea$Rank1 <- 1:nrow(One_drug_Gsea)
  
  after_filter_GSEA_list[[i]] <- One_drug_Gsea
}

names(after_filter_GSEA_list) <- names(mol_result_list_score)
saveRDS(after_filter_GSEA_list, "Step2 data/after_filter_GSEA_list.RDS")

Organ_drug_Gsea_data <- do.call(rbind, after_filter_GSEA_list)


##### Combine Random Walk and GSEA Scores into a Weighted Score #####
mol_result_list_score <- mol_result_list_score[-4]
merged_df_list <- list()

for (i in 1:length(after_filter_GSEA_list)) {
  ### GSEA drug score dataframe
  One_drug_GSEA <- after_filter_GSEA_list[[i]]
  
  ### Random walk drug score dataframe
  One_drug_Random <- mol_result_list_score[[i]]
  One_drug_Random$Rank2 <- 1:nrow(One_drug_Random)
  colnames(One_drug_Random)[1:2] <- c("drug", "Random_score")
  
  ### Merge GSEA and Random walk scores by drug name
  merged_df <- merge(One_drug_GSEA, One_drug_Random, by = "drug")
  
  ### Normalize ranks to [0,1]
  normalize <- function(x) {
    return((x - min(x)) / (max(x) - min(x)))
  }
  
  merged_df$Standard_NES <- 1 / merged_df$Rank1
  merged_df$Standard_Random <- 1 / merged_df$Rank2
  
  merged_df$Standard_NES <- normalize(merged_df$Standard_NES)
  merged_df$Standard_Random <- normalize(merged_df$Standard_Random)
  
  ### Calculate weighted rank: 40% from GSEA NES, 60% from Random walk score
  merged_df$WeightedRank <- 0.4 * (merged_df$Standard_NES) + 0.6 * merged_df$Standard_Random
  
  ### Sort by weighted rank descending
  merged_df <- merged_df[order(merged_df$WeightedRank, decreasing = TRUE), ]
  merged_df$total_rank <- 1:nrow(merged_df)
  
  merged_df_list[[i]] <- merged_df
}

merged_df_data <- do.call(rbind, merged_df_list)
