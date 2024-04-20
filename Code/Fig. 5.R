# Figure 5

## 5a-b
library(ConsensusClusterPlus)
library(RColorBrewer)
dataframe <- as.data.frame(read.csv("yellow.csv",as.is = TRUE, encoding = 'UTF-8'));rownames(dataframe) <- dataframe[,1];dataframe <- dataframe[,-1]
results = ConsensusClusterPlus(as.matrix(dataframe), maxK=10, reps=1000, pItem=0.8, pFeature=1, tmyPal=c('white','#FDB462'), clusterAlg="km", distance="euclidean", seed=12621, plot="pdf", writeTable=TRUE)

# 5c
library(gsva)
library(pheatmap)
# Load data
uni_matrix <- as.data.frame(read.csv("TPM.csv", as.is = TRUE, encoding = 'UTF-8'))
rownames(uni_matrix) <- uni_matrix[,1]
uni_matrix <- uni_matrix[,-1]
gene_set <- read.csv("mmc3.csv")[, 1:2]
list <- split(as.matrix(gene_set)[,1], gene_set[,2])
# Perform ssGSEA analysis
gsva_matrix <- gsva(as.matrix(uni_matrix), list, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
gsva_matrix1 <- t(scale(t(gsva_matrix)))
gsva_matrix1[gsva_matrix1 < -2] <- -2
gsva_matrix1[gsva_matrix1 > 2] <- 2
# Define gene sets
anti_tumor <- c('Activated CD4 T cell', 'Activated CD8 T cell', 'Central memory CD4 T cell', 'Central memory CD8 T cell', 'Effector memeory CD4 T cell', 'Effector memeory CD8 T cell', 'Type 1 T helper cell', 'Type 17 T helper cell', 'Activated dendritic cell', 'CD56bright natural killer cell', 'Natural killer cell', 'Natural killer T cell')
pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'CD56dim natural killer cell', 'Immature dendritic cell', 'Macrophage', 'MDSC', 'Neutrophil', 'Plasmacytoid dendritic cell')
anti <- gsub('^ ','', rownames(gsva_matrix1)) %in% anti_tumor
pro <- gsub('^ ','', rownames(gsva_matrix1)) %in% pro_tumor
non <- !(anti|pro)
gsva_matrix1 <- rbind(gsva_matrix1[anti,], gsva_matrix1[pro,], gsva_matrix1[non,])
# Normalization function
normalization <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
# Normalize data
nor_gsva_matrix1 <- normalization(gsva_matrix1)
# Generate heatmap
pheatmap(mat = nor_gsva_matrix1,
         show_colnames = F,
         cluster_rows = F, cluster_cols = T,
         clustering_method = "ward.D2",
         fontsize = 10, gaps_row = c(12,20),
         filename = 'ssgsea.pdf', width = 8,
         border_color = NA)

## 5d-e
library(ggplot2)
ggplot(data, aes(logFC, -log10(Pvalue))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999") +
  geom_point(aes(size = -log10(Pvalue), color = -log10(Pvalue))) +
  scale_color_gradientn(values = seq(0, 1, 0.2),colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  scale_size_continuous(range = c(1, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.text = element_text(size = 12),
        legend.position = "left",
        legend.justification = c(0, 1)) +
  ylim(-0, 31) +
  guides(col = guide_colourbar(title = "-Log10(Pvalue)"),
         size = "none") +
  geom_text(aes(label = label, color = -log10(Pvalue)), color = c("red"),
            size = 4, vjust = 1.5, hjust = 1) +
  labs(title = "Validation of mRNA sequence in CAF module", size = 15) +
  xlab("Log2FC") +
  ylab("-Log10(Pvalue)")
ggsave("Validation of mRNA sequence in CAF module.pdf", height = 8.5, width = 10)

## 5g 
library(ggpubr)
library(ggrepel)
timer <- read.csv(file = "timer.csv")
pdf(file = "Dysfunction.pdf", width = 7, height = 7) 
ggplot(timer, aes(x = SYNCRIP, y = Dysfunction)) +
  geom_point(color = "#3B7FC7", size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, color = "black", face = "bold"),
        axis.text.y = element_text(size = 10, color = "black", face = "bold")) +
  geom_smooth(method = 'lm', se = TRUE, show.legend = FALSE, linetype = 1, size = 0.6, color = '#333333') +
  stat_cor(method = "spearman", show.legend = FALSE) +
  xlab('SYNCRIP') +
  ylab("Dysfunction")
dev.off()






















