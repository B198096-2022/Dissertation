---
title: "DoubletFinder_pipeline"
author: "B198096"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


```{r, initializing}
library("DoubletFinder")
#https://github.com/chris-mcginnis-ucsf/DoubletFinder
library(Seurat)

sce_list <- list()
files <- list()

range <- 1:9
dir <- "~/Desktop/Bioinformatics/Dissertation/ramachandran_seurat/"
files <- paste0(dir,"ramachandran_fibrosis_cell_",range,"_UMAP.rds")


for (i in range){
  sce_list[i] <- readRDS(file = files[[i]])
} 


```


```{r, expected doublet calc}
exp_d <- list()
for (i in range){
  #Total cells recovered 
  cells <- ncol(sce_list[[i]])
  #10x calcualtes expected doublets as a function of number of recovered cells 
  #Their table is for brackets of 1000 cells, so I am setting the boundries at 500 cells 
  #So that the bracket used is whatever 1000s is closest 
  bracket <- as.integer(cells/1000)
  bracket_up <- cells%%1000
  #Then this is a conditional to assign expected doublet rate   
  #Based on the bracket information 
  if (cells < 500) {
    exp_d[[i]] <- 0.004
    } else {
      exp_d[[i]] <- 0.008 * bracket
      if (bracket_up > 500) {
        exp_d[[i]] <- 0.008 * (bracket + 1)
      }
    }
  }
```

```{r, pK optimization}
pK <- list()
for (i in range){
  #This simulates various pK values for this data to identify the optimal value 
  #I am specifying that 30 PCs were used in the original PCA (I could change 
  #This for the DF analysis, but choosing to use all PCs)
  #And specifying that I used SCTransform 
  sweep.res.list_tf1 <- paramSweep_v3(sce_list[[i]], PCs = 1:30, sct = TRUE)
  #Generates a summary table of the above simultions. 
  #No multiplexed data so GT is false 
  sweep.stats_tf1 <- summarizeSweep(sweep.res.list_tf1, GT = FALSE)
  #Then this calculates BC metrics as a score of success for each pK value 
  bcmvn_tf1 <- find.pK(sweep.stats_tf1)
  #This goes through and pulls the optimal pK value after the simulations 
  pK_temp <- bcmvn_tf1 %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  pK[[i]] <- as.numeric(as.character(pK_temp[[1]]))
}
```


## Including Plots

You can also embed plots, for example:

```{r expected doublet adjustement}
nExp_poi <- list()
nExp_poi.adj <- list()
for (i in range){
  #Adding the annotations of cluster results from the PCA 
  annotations <- sce_list[[i]]@meta.data$ClusteringResults
  #Modelling homotypic doublets based on clustering patterns 
  homotypic.prop <- modelHomotypic(annotations)
  #Predicting the number of expected doublets based on exp_d from above 
  nExp_poi[[i]] <- round(exp_d[[i]]*nrow(sce_list[[i]]@meta.data))  
  #Adjusting the total expected doublet count
  nExp_poi.adj[[i]] <- round(nExp_poi[[i]]*(1-homotypic.prop))
  }
```


```{r, running DF}
df_col <- list()
for (i in range){
  #Running DoubletFinder 
  #I am setting PCs to 30 again 
  #The developers found that results are not significantly impacted by pN 
  #And so I followed their default value of 0.25 
  #pK and nExp are calculated above 
  #reuse.pANN is not being used, and SCTransform was used
  sce_list[[i]] <- doubletFinder_v3(sce_list[[i]], PCs = 1:30, pN = 0.25, pK = pK[[i]], nExp = nExp_poi[[i]], reuse.pANN = FALSE, sct = TRUE)
  #sce_list[[i]] <- doubletFinder_v3(sce_list[[i]], PCs = 1:30, pN = 0.25, pK = pK[[i]], nExp = nExp_poi.adj[[i]], reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
  df_col[[i]] <- paste0("DF.classifications_0.25_",pK[[i]],"_",nExp_poi[[i]])
  DimPlot(sce_list[[i]], reduction = 'umap', group.by = df_col[[i]])
}
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

```{r}
pK <- list()
nExp_poi <- list()
nExp_poi.adj <- list()
df_col <- list()
for (i in range){
  #This simulates various pK values for this data to identify the optimal value 
  #I am specifying that 30 PCs were used in the original PCA (I could change 
  #This for the DF analysis, but choosing to use all PCs)
  #And specifying that I used SCTransform 
  sweep.res.list_tf1 <- paramSweep_v3(sce_list[[i]], PCs = 1:30, sct = TRUE)
  #Generates a summary table of the above simultions. 
  #No multiplexed data so GT is false 
  sweep.stats_tf1 <- summarizeSweep(sweep.res.list_tf1, GT = FALSE)
  #Then this calculates BC metrics as a score of success for each pK value 
  bcmvn_tf1 <- find.pK(sweep.stats_tf1)
  #This goes through and pulls the optimal pK value after the simulations 
  pK_temp <- bcmvn_tf1 %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  pK[[i]] <- as.numeric(as.character(pK_temp[[1]]))
  #Adding the annotations of cluster results from the PCA 
  annotations <- sce_list[[i]]@meta.data$ClusteringResults
  #Modelling homotypic doublets based on clustering patterns 
  homotypic.prop <- modelHomotypic(annotations)
  #Predicting the number of expected doublets based on exp_d from above 
  nExp_poi[[i]] <- round(exp_d[[i]]*nrow(sce_list[[i]]@meta.data))  
  #Adjusting the total expected doublet count
  nExp_poi.adj[[i]] <- round(nExp_poi[[i]]*(1-homotypic.prop))
  #Running DoubletFinder 
  #I am setting PCs to 30 again 
  #The developers found that results are not significantly impacted by pN 
  #And so I followed their default value of 0.25 
  #pK and nExp are calculated above 
  #reuse.pANN is not being used, and SCTransform was used
  sce_list[[i]] <- doubletFinder_v3(sce_list[[i]], PCs = 1:30, pN = 0.25, pK = pK[[i]], nExp = nExp_poi[[i]], reuse.pANN = FALSE, sct = TRUE)
  #sce_list[[i]] <- doubletFinder_v3(sce_list[[i]], PCs = 1:30, pN = 0.25, pK = pK[[i]], nExp = nExp_poi.adj[[i]], reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
  df_col[[i]] <- paste0("DF.classifications_0.25_",pK[[i]],"_",nExp_poi[[i]])
  DimPlot(sce_list[[i]], reduction = 'umap', group.by = df_col[[i]])
}
```


DimPlot(sce_list[[2]], reduction = 'umap', group.by = df_col[[2]])

