---
title: "QC_Filter_test"
author: "B198096"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

## Setup environment and load required packages
```{r}
#create list of packages to load
#install.packages('SoupX')
#install.packages('bigutilsr')
load_pkg <- rlang::quos(Seurat, #scRNA-seq wrapper
scater, #scRNA-seq QC
SoupX, #ambient RNA removal
bigutilsr, #PCA outlier detection
ggplot2, #visualisation
dplyr
)
#load packages
invisible(lapply(lapply(load_pkg, rlang::quo_name),
library,
character.only = TRUE
))
```

## Load in data
```{r}
rf_cells <- list()
sce_list <- list()
range <- 1:9
out_name <- "ramachandran_fibrosis_filtered.RDS"
dir <- "~/Desktop/Bioinformatics/Dissertation/ramachandran_seurat/"
fig_path <- "~/Desktop/Bioinformatics/Dissertation/ramachandran_seurat/"
files <- paste0(dir,"ramachandran_fibrosis_cell_",range,"_UMAP.rds")
for (i in range){
  sce_list[i] <- readRDS(file = files[[i]])
} 

sce_list[[1]] <- rf1
sce_list[[2]] <- rf2
sce_list[[3]] <- rf3
sce_list[[4]] <- rf4
sce_list[[5]] <- rf5
sce_list[[6]] <- rf6
sce_list[[7]] <- rf7
sce_list[[8]] <- rf8
sce_list[[9]] <- rf9

#object_merged <- merge(rf1, y = c(rf2, rf3, rf4, rf5, rf6, rf7, rf8, rf9), add.cell.ids = c("rf1", "rf2", "rf3", "rf4", "rf5", "rf6", "rf7", "rf8", "rf9"), project = "rf_merged_QC")

#Save the unfiltered dimensions for expectd doublet values
for (i in range){
  rf_cells[[i]] <- ncol(sce_list[[i]])
}
saveRDS(rf_cells, file = paste0(dir,"rf_unfiltered_dim"))

```


## Filtered counts matrix for read outliers
```{r}
#object_merged <- merge(rf1, y = rf2, add.cell.ids = c("rf1", "rf2"), project = "rf_big")


#determine QC outliers
#https://www.rdocumentation.org/packages/scater/versions/1.0.4/topics/isOutlier
sce <- as.SingleCellExperiment(object_merged)
sce <- addPerCellQC(sce, subsets=list(mito=grep("^MT-", rownames(sce))))
count_outlier <- isOutlier(sce$detected, nmads=3, type="both", log=TRUE, batch = sce$broad_enrichment)
detected_outlier <- isOutlier(sce$sum, nmads=3, type="both", log=TRUE, batch = sce$broad_enrichment)
mt_outlier <- isOutlier(sce$subsets_mito_percent, nmads=5, type="higher", batch = sce$broad_enrichment)
sce$count_outlier <- count_outlier
sce$detected_outlier <- detected_outlier
sce$mt_outlier <- mt_outlier
sce$final_outlier <- count_outlier | detected_outlier | mt_outlier
#plot QC outliers
plt <- plotColData(sce, x = "sum", y="detected", colour_by="final_outlier")
ggsave(filename = paste0(fig_path, "sum_detected_scatter_by_counts_outlier.pdf"), plot = plt)
plt
plt <- plotColData(sce, x = "subsets_mito_percent", y="detected", colour_by="final_outlier")
ggsave(filename = paste0(fig_path, "pct_mt_detected_scatter_by_counts_outlier.pdf"), plot = plt)
plt
#remove QC outliers
dim(sce)
sce2 <- sce[,which(!sce$final_outlier)]
dim(sce2)
object_merged2 <- subset(object_merged, cells = colnames(sce2))
```


## Filter counts matrix for expression outliers
```{r}
#determine expression outliers
object_merged2 <- NormalizeData(object_merged2)
VariableFeatures(object_merged2) <- SelectIntegrationFeatures(lapply(X = SplitObject(object_merged2, split.by = "orig.ident"), FUN = function(x) {x <- FindVariableFeatures(x, verbose = F)}))
object_merged2 <- object_merged2 %>%
ScaleData() %>%
RunPCA()
x <- as.matrix(object_merged2[["pca"]]@cell.embeddings)
y <- LOF(x, log = F, ncores = 8)
outliers <- isOutlier(y, log = T, type = "higher")
object_merged2$expression_outliers <- outliers
#remove expression outliers
dim(object_merged2)
object_merged3 <- subset(object_merged2, expression_outliers == FALSE)
dim(object_merged3)
#rm(object_merged_list)
```


## Filtered counts matrix for weakly expressed genes
```{r}
#remove weakly expressed genes
counts <- GetAssayData(object = object_merged3, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) > 3
removed_genes <- setdiff(rownames(object_merged3), names(keep_genes))
filtered_counts <- counts[keep_genes, ]
dim(object_merged3)
#clean object meta data
#Create a new final object with these filtered counts from above 
object_merged4 <- CreateSeuratObject(filtered_counts) 
                                     
meta.data = object_merged3@meta.data[,c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_decontXcounts", "nFeature_decontXcounts", "expression_outliers", "log10GenesPerUMI", "percent.rb")])
dim(object_merged4)
object_merged4$log10GenesPerUMI <- log10(object_merged4$nFeature_RNA) / log10(object_merged4$nCount_RNA)
object_merged4$percent.mt <- PercentageFeatureSet(object = object_merged4, pattern = "^MT-")
object_merged4$percent.rb <- PercentageFeatureSet(object = object_merged4, pattern = "^RP[SL]")
object_merged4@meta.data[, colnames(meta_data)] <- ""
for(i in rownames(meta_data)){
object_merged4@meta.data[which(object_merged4@meta.data$orig.ident == i), colnames(meta_data)] <- meta_data[i,]
}
dim(object_merged4)
saveRDS(object_merged4, file = paste0(dir, "outname.RDS"))
```

"orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_decontXcounts", "nFeature_decontXcounts", "expression_outliers", "log10GenesPerUMI", "percent.rb" 

