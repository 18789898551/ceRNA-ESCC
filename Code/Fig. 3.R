# Figure 3

library(WGCNA)
library(edgeR)
library(tidyverse)
gene<-read.csv("TPM.csv",row.names = 1,header = T)
datExpr0 = as.data.frame(gene) # 行是基因，列是样本
# Load clinical information
allTraits = read.csv("clinical.csv",row.names = 1,header = T)
# Match expression level with trait data
fpkmSamples = rownames(gene) 
traitSamples =rownames(allTraits) 
traitRows = match(fpkmSamples, traitSamples) 
datTraits = allTraits[traitRows,] 
collectGarbage()

## 3b
pdf(file = "2.Sample dendrogram and trait heatmap.pdf", width = 15, height = 8) 
sampleTree2 = hclust(dist(gene), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) #用颜色代表关联度
plotDendroAndColors(sampleTree2, traitColors,
                    cex.colorLabels = 0.9, 
                    cex.dendroLabels = 0.6, 
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Soft threshold
powers = c(c(1:10), seq(from = 10, to=15, by=2))
# Obtain R square and average connectivity under each threshold
sft = pickSoftThreshold(gene, powerVector = powers, verbose = 5)
sft$powerEstimate

datExpr0[] <- lapply(datExpr0, as.numeric)
str(datExpr0)
datExpr0 <- as.data.frame(datExpr0)
cor <- WGCNA::cor
net = blockwiseModules(datExpr0, power = 5,
                       TOMType = "unsigned", 
                       minModuleSize = 30,
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
cor<-stats::cor

# Save the allocation module and the gene information contained in the module
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, file = "networkConstruction.RData")
### Module correlates phenotypic data and identifies important genes
### Module-phenotype data association
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
### Recalculate modules with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
### Color-code each association by its associated value
### Display the correlation coefficient and P value between the module and the phenotypic data
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

## 3e
pdf(file = "5.Module-trait relationships.pdf", width = 12, height = 6.4) 
par(mar = c(7, 8.5, 5, 5))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               font.lab.y = 2,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

## 3c
pdf(file = "7.Eigengene adjacency heatmap.pdf", width =5, height = 4) 
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                      plotHeatmaps = TRUE, 
                      plotDendrograms = T,
                      colorLabels = TRUE, 
                      
                      plotAdjacency = F,
                      greyLabel = "grey", 
                      marHeatmap = c(3,4,2,2), 
                      xLabelsAngle = 90)

dev.off()

library(corrplot)
library(pheatmap)

## 3d
pdf(file="stemness.pdf", width=3, height=3)
stemness <- read.csv("correlation_stemness.csv",row.names = 1,header = T)#样本
ciber.res<-stemness[,colSums(stemness) > 0]
cor_ciber<- cor(ciber.res, method = 'pearson') 
cor_ciber <- round(cor(ciber.res), 2)
col= colorRampPalette(c( '#2B8CBE', 'white','#EF3B2C' ))(40)
res1 <-cor.mtest( ciber.res, conf.level= .95)#显著水平
corrplot(cor_ciber,method = "color",
         title = "",tl.col = "black",tl.cex = 0.7,rect.col = "black",
         type= "upper",col = col,p.mat = res1$p,
         sig.level = c(.001, .01, .05),insig = "label_sig",pch.cex = 0.5)
dev.off()
