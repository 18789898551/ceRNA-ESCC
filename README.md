# ceRNA-ESCC
Immunotherapy for esophageal squamous cell carcinoma (ESCC) exhibits notable variability in efficacy. Concurrently, recent research emphasizes circRNAs' impact on the ESCC tumor microenvironment. To further explore the relationship, we leveraged circRNA, microRNA, and mRNA sequence datasets to construct a comprehensive immune-related circRNA–microRNA–mRNA network, revealing competing endogenous RNA (ceRNA) roles in ESCC. The network comprises 16 circular RNAs, 13 microRNAs, and 1,560 mRNAs. Weighted gene co-expression analysis identified immune-related modules, notably cancer-associated fibroblast (CAF) and myeloid-derived suppressor cell modules, correlating significantly with immune and stemness scores. Among them, the CAF module plays a crucial role in extracellular matrix function and effectively discriminates ESCC patients. Four hub collagen family genes within CAF correlated robustly with CAF, macrophage infiltration, and T-cell exclusion. In-house sequencing and RT-qPCR validated their elevated expression. We also identified CAF module-targeting drugs as potential ESCC treatments. In summary, we established an immune-related circRNA–miRNA–mRNA network that not only illuminates ceRNA functionality but also highlights circRNAs' involvement in the CAF through collagen gene targeting. These findings hold promise to predict ESCC immune landscapes and therapy responses, ultimately aiding in more personalized and effective clinical decision-making.

![Flow chart](https://github.com/18789898551/ceRNA-ESCC/blob/main/Code/Flowchart.png)

## Fig. 2.R: Identification of DECs, DEMs, and DEGs in GEO and TCGA datasets.
### 2a-b:Identification of robust different expressed genes in GEO datasets through RRA.
The analysis identified differentially expressed genes across three distinct datasets, pinpointing commonly upregulated and downregulated genes. Initially, the necessary libraries and data files were imported. Subsequently, each dataset was processed to extract lists of upregulated and downregulated genes. The RobustRankAggreg package was then employed to rank and aggregate these genes, identifying those significantly upregulated or downregulated. Finally, these genes were consolidated into a single table, and a heatmap was generated to visualize their expression patterns across the different datasets.
### 2c: Heatmap depicts the expression levels of the top 100 DEGs in TCGA-ESCC dataset
For heatmap generation to display gene expression data, the process began with the importation of essential R packages and datasets, followed by the establishment of a color palette. Subsequent preprocessing steps included replacing zero values with 0.01 times the minimum value and eliminating rows with a standard deviation of zero. Finally, the pheatmap function was utilized to create the heatmap, which was then saved as a PDF file.

## Fig. 3.R: Characteristic mRNAs selection.
### 3b: Sample dendrogram and trait heatmap in TGCA-ESCC dataset.
The code performs hierarchical clustering of samples using the average linkage method and assigns colors to each sample based on clinical trait data to indicate correlation. It then generates a combined plot of the sample dendrogram and trait heatmap.
### 3e：Module-trait relationship.
Calculates the correlation coefficients and P-values between the eigengenes of each module and the clinical trait data, and presents the results as a heatmap where rows represent modules, columns represent clinical traits, and the color intensity indicates the magnitude of the correlation coefficients, with corresponding P-values displayed alongside.
### 3c: Eigengene adjacency heatmap was depicted to depict the connectivity between different modules in the dataset.
Eigengene adjacency heatmap, displaying the correlations among module eigengenes, and incorporating a hierarchical clustering dendrogram of the modules.
### 3d: Chorograms were derived to illustrate the relationship between immune infiltration and stemness in the dataset.
The Pearson correlation coefficients between the stemness index and other variables were calculated, a multiple test was performed on the correlation matrix to determine the significance level, and a correlation matrix plot was generated, incorporating color-coded correlation coefficients and significance labels.

## Fig. 4.R: ceRNA construction and immune-related module exploration.
## 4b/d: The GO function enrichment of DEGs.
A Gene Ontology (GO) enrichment analysis was performed on a designated gene set to identify the foremost 10 terms in each of the three ontology categories: biological process (BP), cellular component (CC), and molecular function (MF). The results were subsequently aggregated into a structured data frame. Following this, a bar chart depicting GO enrichment was created utilizing the ggplot2 package, offering a clear graphical demonstration of the gene count distributions among the varied ontological domains.

## Fig. 5.R: Investigating Immune Infiltration in the Key Cluster of the CAF Module and Validating Hub Gene Expression.
### 5a-b: The CDF curves of consensus matrix for each k (indicated by colors) and consensus matrix heatmaps
The ConsensusClusterPlus package was utilized to perform consensus clustering analysis on the data, with the maximum number of clusters set to 10 and the number of repetitions set to 1000. The k-means algorithm and Euclidean distance metric were employed, with various parameters configured and the output type specified as a PDF file. Finally, a data table of the clustering results was generated.
### 5c: Comparison of the relationships between the infiltration abundance of 28 immune cell subsets of two clusters. 
The gsva package was employed to perform single-sample gene set enrichment analysis (ssGSEA), followed by the utilization of the pheatmap package to generate a heatmap illustrating the enrichment profiles of various gene sets across samples. The data was then subjected to normalization and clustering analysis.
### 5d-e: Volcano maps separately showing DECs and DEGs between tumor and normal samples
The ggplot2 package was utilized to create a scatter plot that displays the validation results of mRNA sequences in the Cancer-Associated Fibroblast module. The size and color of the points were adjusted to represent statistical significance. Relevant legends, titles, and axis labels were added to the plot, which was subsequently saved as a PDF file.
### 5g: Scatterplots between 4 hub genes expression in CAF module with Macrophage. 
A scatter plot was employed to illustrate the relationship between SYNCRIP and Dysfunction. The plot incorporated a linear regression model for fitting and included the Spearman correlation coefficient as a statistical measure. Details such as plot titles, axis labels, and font sizes were adjusted accordingly.

## Fig. 6.R: Immune correlation of hub genes in MDSC module and potential compound identified in CAF module.
### 6a-b: Scatterplots between 4 hub genes expression in CAF module
A scatter plot that demonstrates the relationship between SYNCRIP and Dysfunction. 
### 6c: Exploration of the distribution of four hub genes using single-cell sequencing.
(Left)The single-cell RNA sequencing data was preprocessed, subjected to dimensionality reduction, and clustering analysis. Subsequently, the UMAP algorithm was employed for visualization. Finally, cell types were annotated using the HumanPrimaryCellAtlasData and SingleR packages.

(right)The FeaturePlot function was utilized to visualize the expression of specific genes (COL1A1, COL5A1, COL5A2, and COL4A1) in a single-cell RNA sequencing dataset, applying a designated color coding and border theme. 

### We have uploaded the code corresponding to each figure image to facilitate future study and discussion.
