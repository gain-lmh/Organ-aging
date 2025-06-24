# Load required libraries
library(GSEABase)
library(GSVA)
library(ggplot2)
library(cowplot)
library(ggridges)
library(Hmisc)
library(edgeR)
library(magrittr)
library(gtools)
library(RColorBrewer)

# === Configuration ===
setwd("/home/data/t140311/Organ_aging")
dir_list <- dir()

# Load gene sets
age_genes <- read.table("../data/share_data/group_time_gene.txt", header = TRUE, sep = "\t")
group_organ_num <- read.table("/home/data/t140311/Organ_aging/share_data/group_gene_Organ_num.txt", header = TRUE)
up_genes <- age_genes[age_genes$regulation == "up_gene", ]
down_genes <- age_genes[age_genes$regulation == "down_gene", ]

# Storage lists
ridge_plots <- list()
cor_plots <- list()
bar_plots <- list()
sample_scores <- list()
group_scores <- list()

# Normalize score vector
normalize_score <- function(score, up_n, down_n) {
  normalized <- score / (sd(score) * sqrt((1 / up_n) + (1 / down_n)))
  return(normalized / 5)
}

# Loop through organs
total_organs <- length(dir_list)
for (j in seq_len(total_organs)) {
  setwd(paste0("D:/AgingStudy/Result/", dir_list[j]))
  
  cluster_data <- readRDS("consensus_cluster.rds")
  exp_data <- log2(cpm(as.matrix(readRDS("all_data_exp.rds"))) + 1)
  exp_data <- as.data.frame(t(exp_data))

  # Calculate up-regulated gene score
  up_exp <- exp_data[, up_genes$gene, drop = FALSE]
  up_weight <- group_organ_num[group_organ_num$gene %in% up_genes$gene, ]$Organ_num * 0.1
  up_exp <- sweep(up_exp, 2, up_weight, `*`)
  up_score <- rowMeans(up_exp)

  # Calculate down-regulated gene score
  down_exp <- exp_data[, down_genes$gene, drop = FALSE]
  down_weight <- group_organ_num[group_organ_num$gene %in% down_genes$gene, ]$Organ_num * 0.1
  down_exp <- sweep(down_exp, 2, down_weight, `*`)
  down_score <- rowMeans(down_exp)

  final_score <- normalize_score(up_score - down_score, ncol(up_exp), ncol(down_exp))

  cluster_data$Time <- sapply(strsplit(as.character(cluster_data$Time), split = "_"), `[`, 1)
  cluster_data$Time <- as.numeric(mixedsort(cluster_data$Time))

  df <- data.frame(Time = cluster_data$Time, Score = final_score, Point = as.factor(cluster_data$Time))

  # Ridge plot
  colors <- c("#9CE1BB", "#8DEA86", "#B76BCC", "#CCD880", "#A8D9E3", "#E4BDD1",
              "#E8957F", "#9AA5DA", "#E086C3", "#E0DCC1")
  ridge_plot <- ggplot(df, aes(x = Score, y = Point, fill = Point, color = Point)) +
    theme_classic() +
    geom_density_ridges(alpha = 1, scale = 2, show.legend = FALSE) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    ylab(dir_list[j]) + xlab("ssGSEA score") +
    theme(axis.text.y = element_text(face = "italic"), axis.text.x = element_text(face = "bold"))
  ridge_plots[[j]] <- ridge_plot

  # Correlation plot
  cor_res <- rcorr(df$Time, df$Score)
  cor_val <- round(cor_res$r[1, 2], 2)
  p_val <- ifelse(cor_res$P[1, 2] > 0.01, formatC(cor_res$P[1, 2], digits = 2), "p < 0.01")
  cor_plot <- ggplot(df, aes(Time, Score)) +
    geom_point(color = "#988d7b") +
    geom_smooth(method = "lm", fill = "#b2e7fa", color = "#00aeef", alpha = 0.8) +
    scale_x_continuous(breaks = c(1, seq(3, 27, 3))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(face = "bold.italic"),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = paste0(dir_list[j], "  ρ =", cor_val, "  ", p_val))
  cor_plots[[j]] <- cor_plot

  # Bar plot
  bar_data <- aggregate(Score ~ Time, df, mean)
  bar_data$Time <- 1:10
  organ_colors <- c("#866DA6", "#38A8D7", "#E8C1C7", "#CD3E39", "#EEC684", "#BAD09E", "#D77DA7",
                    "#CEC1D7", "#A7CCE5", "#B6773D", "#8EB69C", "#98D7EA", "#E1D572",
                    "#E7AFC9", "#67B1C9", "#189A52")
  bar_plot <- ggplot(bar_data, aes(Time, Score)) +
    geom_bar(stat = "identity", fill = organ_colors[j], color = "gray2", width = 0.75) +
    geom_point(color = "blue", size = 4) +
    geom_line(color = "blue", size = 1.3) +
    scale_x_continuous(breaks = 1:10, labels = c(1, seq(3, 27, 3))) +
    theme_classic() +
    ggtitle(dir_list[j]) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(face = "bold.italic"))
  bar_plots[[j]] <- bar_plot

  # Save scores
  df$Group <- cluster_data$group
  group_avg <- aggregate(Score ~ Group, df, mean)
  group_avg$Organ <- paste(dir_list[j], group_avg$Group, sep = "_")
  group_scores[[j]] <- group_avg
  sample_scores[[j]] <- df
}

# Save results
saveRDS(sample_scores, "../../data/GTAG_sample_score_list.RDS")

# Combine plots
plot_grid_all <- function(plot_list, output_file) {
  combined <- plot_grid(plotlist = plot_list, ncol = 4)
  ggsave(filename = output_file, plot = combined, height = 10, width = 14)
}

plot_grid_all(ridge_plots, "/home/data/t140311/Organ_aging/Figures/ridge_all_organ.pdf")
plot_grid_all(cor_plots, "/home/data/t140311/Organ_aging/Figures/correlation_all_organ.pdf")
plot_grid_all(bar_plots, "/home/data/t140311/Organ_aging/Figures/bar_all_organ.pdf")

# Polar plot
summary_data <- do.call(rbind, group_scores)
summary_data$x <- summary_data$Score + 26
summary_data$label <- summary_data$x - 32

ggplot(summary_data, aes(x = Organ, y = x, fill = Group)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = label), color = "white", vjust = 1.5) +
  coord_polar(theta = "x") +
  scale_fill_manual(values = c("1" = "#809EBA", "2" = "#F8BBBB", "3" = "#AC9B7B")) +
  ylim(c(-1, 46)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, color = "black"))

ggsave("/home/data/t140311/Organ_aging/Figures/group_gene_gsea_polar.pdf", width = 9, height = 8)
