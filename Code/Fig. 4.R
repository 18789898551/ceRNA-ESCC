# Figure 4

## 4b/d
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
library(clusterProfiler)
library(GOplot)
library(DOSE)
library(ggnewscale)
library(topGO)
library(circlize)
library(ComplexHeatmap)
### GO enrichment
### CC: cellular components; MF: molecular function; BP: biological process
ego_CC <- enrichGO(gene = id,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = TRUE)
ego_BP <- enrichGO(gene = id,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = TRUE)
ego_MF <- enrichGO(gene = id,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = TRUE)
### Filter the top 10 processes
display_number = c(10, 10, 10)
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

### Reassemble into data frame
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID), 
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))
for(i in 1:nrow(go_enrich_df)){
  description_splite=strsplit(go_enrich_df$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
  go_enrich_df$Description[i]=description_collapse
  go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
}

## Draw GO histogram
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#F38764","#F8E639","#D65346")
library(ggplot2)
pdf(file = "GO.pdf", width = 10, height = 7) 
ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) +
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=13))
dev.off()




























