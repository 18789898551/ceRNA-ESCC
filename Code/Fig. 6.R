# Figure 6

## 6a-b
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

## 6c
rm(list = ls())
load("D:\\single cell\\GSE196756\\harmony.RData")
### Visualize content of RNA-gene 
plot1 <- FeatureScatter(scRNA_harmony, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA_harmony, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
m_fit <- subset(scRNA_harmony, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25)
m_fit_norm <- NormalizeData(m_fit , normalization.method = "LogNormalize", scale.factor = 10000)
### Scale data
all.genes <- rownames(m_fit_norm) 
m_fit_norm<-FindVariableFeatures(m_fit_norm) #找到所有特征（可不执行）
m_fit_norm<- ScaleData(m_fit_norm)
m_fit_norm_sca<- ScaleData(m_fit_norm, features = all.genes)
# PCA dimensionality reduction
m_fit_norm_pca <- RunPCA(m_fit_norm_sca,
                         features = VariableFeatures(object = m_fit_norm_sca))
### UMAP Visualization
### Find the best number of clusters first
m_fit_norm_pca_c <- FindNeighbors(m_fit_norm_pca, dims = 1:10)
m_fit_norm_pca_c <- FindClusters(m_fit_norm_pca_c, resolution = 0.5)
UMAP <- RunUMAP(m_fit_norm_pca_c, dims = 1:10)
DimPlot(UMAP, reduction = "umap")
### Annotation
meta=UMAP@meta.data 
hpca.se <- HumanPrimaryCellAtlasData()
pbmc_for_SingleR <- GetAssayData(UMAP, slot="data")
pbmc.hesc <- SingleR(test = pbmc_for_SingleR, ref = hpca.se, labels = hpca.se$label.main, 
                     clusters = UMAP@active.ident, de.method="wilcox") 
new.cluster.ids <- pbmc.hesc$pruned.labels
names(new.cluster.ids) <- levels(UMAP)
testname <- RenameIdents(UMAP, new.cluster.ids)
### Annotation results
col <- c("#FABB6E", "#FC8002", "#ADDB88", "#369F2D", "#FAC7B3","#EE4431","#B9181A",
         "#CEDFEF", "#92C2DD", "#4995C6", "#1663A9", "#B4B4D5","#8481BA","#614099")
Figure7 <- DimPlot(testname,label = T, pt.size = 1, cols= col , label.size = 4)+
  NoLegend()+
  labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.ticks.x = element_blank())

pdf(file = "单细胞注释结果.pdf",width =8,height = 7)
plot_grid(Figure7)
dev.off()

# 4 hub genes
color <- c('lightgrey', '#FABB6E','#FC8002')  
Figure3 <- FeaturePlot(UMAP, features = 'COL1A1',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))#加边框
Figure4 <- FeaturePlot(UMAP, features = 'COL5A1',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))#加边框
Figure5 <- FeaturePlot(UMAP, features = 'COL5A2',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))#加边框
Figure6 <- FeaturePlot(UMAP, features = 'COL4A1',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))#加边框

pdf(file = "Gene.pdf",width =12,height = 10)
plot_grid(Figure3,Figure4,Figure5,Figure6)
dev.off()















