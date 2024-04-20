# Figure 2

## 2a-b
rm(list = ls())
options(stringsAsFactors = F)
library(RobustRankAggreg)
library(clusterProfiler)
### Differenced gene list
padj=0.05
logFC=1
### Importing files
files = c("GSE131969.csv","GSE150476.csv","GSE103104.csv")# circRNA-DECs
files = c("GSE114110.csv","GSE43732.csv", "GSE55856.csv" )# miRNA-DEMs
upList=list()
downList=list()
allFCList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt= read.csv(inputFile)
  header=unlist(strsplit(inputFile,"_"))
  downList[[header[1]]]=as.vector(rt[,1])
  upList[[header[1]]]=rev(as.vector(rt[,1]))
  fcCol=rt[,1:2]
  colnames(fcCol)=c("Gene",header[[1]])
  allFCList[[header[1]]]=fcCol
}
mergeLe=function(x,y){
  merge(x,y,by="Gene",all=T)}
newTab=Reduce(mergeLe,allFCList)
rownames(newTab)=newTab[,1]
newTab=newTab[,2:ncol(newTab)]
newTab[is.na(newTab)]=0
### Screen commonly upregulated genes
library(RobustRankAggreg)
upMatrix = rankMatrix(upList)
upAR = aggregateRanks(rmat=upMatrix)
colnames(upAR)=c("Name","Pvalue")
upAdj=p.adjust(upAR$Pvalue,method="BH")
upXls=cbind(upAR,adjPvalue=upAdj)
upFC=newTab[as.vector(upXls[,1]),]
upXls=cbind(upXls,logFC=rowMeans(upFC))
write.table(upXls,file="up_.xls",sep="\t",quote=F,row.names=F)
upSig=upXls[(upXls$Pvalue<padj & upXls$logFC>logFC),]
write.table(upSig,file="upSig_.xls",sep="\t",quote=F,row.names=F)
### Screen commonly down-regulated genes
downMatrix = rankMatrix(downList)
downAR = aggregateRanks(rmat=downMatrix)
colnames(downAR)=c("Name","Pvalue")
downAdj=p.adjust(downAR$Pvalue,method="bonferroni")
downXls=cbind(downAR,adjPvalue=downAdj)
downFC=newTab[as.vector(downXls[,1]),]
downXls=cbind(downXls,logFC=rowMeans(downFC))
write.table(downXls,file="down_.xls",sep="\t",quote=F,row.names=F)
downSig=downXls[(downXls$Pvalue<padj & downXls$logFC< -logFC),]
write.table(downSig,file="downSig_.xls",sep="\t",quote=F,row.names=F)
### Merge commonly up- and down-regulated genes
allSig = rbind(upSig,downSig)
colnames(allSig)
allSig = allSig[,c("Name","logFC")]
write.table(allSig,file = '3个数据集.xls',sep = '\t',quote = F)
### logFC.pdf
hminput=newTab[c(as.vector(upSig[1:20,1]),as.vector(downSig[1:10,1])),]
colnames(hminput) = c("GSE131969","GSE150476","GSE103104")# circRNA
colnames(hminput) = c("GSE114110","GSE43732", "GSE55856" )# miRNA
library(pheatmap)
pdf(file="circRNA_RRA.pdf",width = 8,height = 10)
pheatmap(hminput,display_numbers = TRUE,
         name = "RRA_circRNA",
         fontsize_row=12,
         fontsize_col=12,
         cellwidth = 100,
         cellheight = 15,
         color = colorRamp2(c(-4,-2,0,2,4), c("#0868AC","#4EB3D3","#EEEEEE","#EF3B2C","#CB181D")),
         cluster_cols = FALSE,cluster_rows = FALSE)
dev.off()

## 2c
pkgs <- c("matrixStats", "pheatmap", "RColorBrewer", "tidyverse","ComplexHeatmap","cowplot","ggpubr","bslib","ggthemes","gplots")
lapply(pkgs, library, character.only = T)
### Import samples
seq <- read.csv("gene_.csv",row.names = 1,header = T)#样本
colname <- read.csv("colname.csv",row.names = 1,header = T)#列释义
### Set color
col = list(Type = c("Normal" = "#377EB8", "Tumor" = "#E41A1C"))
seq[seq==0] <- NA
seq[is.na(seq)] <- min(seq,na.rm = T)*0.01
seq <- seq[apply(seq,1, function(x) sd(x)!=0),]
pheatmap(log10(seq))
### Plot
rownames(colname) <- colnames(seq)
pdf(file="heatmap.pdf", width=15, height=15)
pheatmap(as.matrix(log10(seq+1)),
         scale = "row",border_color = "grey70",border=FALSE,
         name = "expression",cellwidth = 6,
         show_colnames = F, show_rownames = T, 
         treeheight_row=20,treeheight_col=20,
         annotation_col  = colname,clustering_method = "ward.D",
         row_km = 2,fontsize = 18,fontsize_row = 8,
         legend = T,cluster_rows = T,cluster_cols = T)
dev.off()
















