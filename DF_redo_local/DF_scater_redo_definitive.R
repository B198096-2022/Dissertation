# DoubletFinder Eddie

#qlogin -l h_vmem=8g
#module load anaconda/5.3.1
#source activate max.Seurat4

#These two packages are used to download old package versions 
library("remotes")
library("fs")

#.libPaths() tells you where your packages are stored 
#"/Library/Frameworks/R.framework/Versions/4.2/Resources/library"
normal_lib <- "/Library/Frameworks/R.framework/Versions/4.2/Resources/library"
#I created a new directory that branches from the standard package directory
old_lib <- "/Library/Frameworks/R.framework/Versions/4.2/Resources/old-versions"

#Then I used the below commands to install the old Seurat and SeuratObject 
#install_version("Seurat", version = "4.3", lib = old_lib)
#install_version("SeuratObject", version = "4.1.3", lib = old_lib)
library("SeuratObject", lib.loc = old_lib)
library("Seurat", lib.loc = old_lib)


options(future.globals.maxSize = 1e12)
library(dplyr)
library(patchwork)
library(DoubletFinder)
#library(tidyverse)
library(ggplot2)

setwd(setwd("/Volumes/cmvm/scs/groups/pramacha-GROUP/Max_Hammer"))

range <- 1:2
sample_author <- "tamburini"
sample_name <- "tamburini_fibrosis_"
sample_label <- "tfc"
fig_lab <- sample_label

datadir <- paste0("/Volumes/cmvm/scs/groups/pramacha-GROUP/Max_Hammer/scater/")
outdir <- "~/Desktop/big_pilot/tamburini/"
#counts <- readRDS(paste0(datadir,sample_author,"/",sample_label,"_raw_counts.rds"))

mt_cutoff <- 30
nfeature_cutoff <- 100
ncount_cutoff <- 0

Files <- list()
Destination <- list()
for (i in range){
  Files[[i]] <- paste0(datadir,sample_author,"/",sample_label, i,"_scater_object.rds")
  Destination[[i]] <- paste0(outdir,sample_name,i,"_DF.rds")
}

Sce <- list()
counts <- list()
for (i in range){
  Sce[[i]] <- readRDS(Files[[i]])
  #In case this is being re-run, Setting default assay back to RNA
  DefaultAssay(Sce[[i]]) <- "RNA"
  #Run SCTransform  
  #Sce[[i]] <- SCTransform(Sce[[i]], verbose = FALSE)
  counts[i] <- dim(Sce[[i]])[2]
  Sce[[i]] <- subset(Sce[[i]], subset = percent.mt < mt_cutoff)
  Sce[[i]] <- subset(Sce[[i]], subset = nFeature_RNA > nfeature_cutoff)
  Sce[[i]] <- subset(Sce[[i]], subset = nCount_RNA > ncount_cutoff)
}

for (i in range){
  Sce[[i]] <- NormalizeData(Sce[[i]])
  Sce[[i]] <- FindVariableFeatures(Sce[[i]], selection.method = "vst", nfeatures = 2000)
  Sce[[i]] <- ScaleData(Sce[[i]])
  #Run PCA
  Sce[[i]] <- RunPCA(Sce[[i]])
  #Run UMAP
  Sce[[i]] <- RunUMAP(Sce[[i]], dims = 1:30)
}


# Run DF
exp_d <- list()
for (i in range){
  #Total cells recovered 
  cells <- as.integer(counts[i])
  #cells <- 130
  #10x calcualtes expected doublets as a function of number of recovered cells 
  #Their table is for brackets of 1000 cells, so I am setting the boundaries at 500 cells
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

#Then Run DoubletFinder, which will label expected doublets  

pK <- list()
nExp_poi <- list()
nExp_poi.adj <- list()
df_col <- list()
for (i in range ){
  #This simulates various pK values for this data to identify the optimal value 
  #I am specifying that 30 PCs were used in the original PCA (I could change 
  #This for the DF analysis, but choosing to use all PCs)
  #And specifying that I used SCTransform 
  sweep.res.list_tf1 <- paramSweep_v3(Sce[[i]], PCs = 1:30, sct = FALSE)
  saveRDS(sweep.res.list_tf1, file = paste0(outdir,fig_lab,"_DF_sweep_res_list.rds"))
  #Generates a summary table of the above simultions. 
  #No multiplexed data so GT is false 
  sweep.stats_tf1 <- summarizeSweep(sweep.res.list_tf1, GT = FALSE)
  saveRDS(sweep.stats_tf1, file = paste0(outdir,fig_lab,"_DF_sweep_stats.rds"))
  #Then this calculates BC metrics as a score of success for each pK value 
  bcvals <- c()
  bcstats <- c()
  for (j in 1:23){
    bcvals <- c()
    for (n in 1:6){
      position <- (j + 23*(n -1))
      #Sometimes DF returns NA for one of the values
      #this removes any potential NA values 
      if (!is.na(sweep.stats_tf1[position, 3])){
        bcvals <- c(bcvals, sweep.stats_tf1[position, 3])
      }
    }
    bcmean <- mean(bcvals)
    bcvar <- var(bcvals)
    bcstat <- bcmean/bcvar
    bcstats <- c(bcstats, bcstat)
    if (bcstat >= max(bcstats)){
      max_pK_pos <- j
    }
  }
  #This goes through and pulls the optimal pK value after the simulations 
  pK_temp <- sweep.stats_tf1[max_pK_pos, 2]
  pK[[i]] <- as.numeric(as.character(pK_temp[[1]]))
  #Adding the annotations of cluster results from the PCA 
  annotations <- Sce[[i]]@meta.data$ClusteringResults
  #Modelling homotypic doublets based on clustering patterns 
  homotypic.prop <- modelHomotypic(annotations)
  #Predicting the number of expected doublets based on exp_d from above 
  nExp_poi[[i]] <- round(exp_d[[i]]*nrow(Sce[[i]]@meta.data))  
  #Adjusting the total expected doublet count
  nExp_poi.adj[[i]] <- round(nExp_poi[[i]]*(1-homotypic.prop))
  #Running DoubletFinder 
  #I am setting PCs to 30 again 
  #The developers found that results are not significantly impacted by pN 
  #And so I followed their default value of 0.25 
  #pK and nExp are calculated above 
  #reuse.pANN is not being used, and SCTransform was used
  Sce[[i]] <- doubletFinder_v3(Sce[[i]], PCs = 1:30, pN = 0.25, pK = pK[[i]], nExp = nExp_poi[[i]], reuse.pANN = FALSE, sct = FALSE)
  #Sce[[i]] <- doubletFinder_v3(Sce[[i]], PCs = 1:30, pN = 0.25, pK = pK[[i]], nExp = nExp_poi.adj[[i]], reuse.pANN = paste0("pANN_0.25_",pK[[i]],"_",nExp_poi[[i]]), sct = FALSE)
  df_col[[i]] <- paste0("DF.classifications_0.25_",pK[[i]],"_",nExp_poi[[i]])
  #DimPlot(Sce[[i]], reduction = 'umap', group.by = df_col[[i]])
}

saveRDS(pK, file = paste0(outdir,fig_lab,"_DF_pK.rds"))

df_results <- list()
#Rename the doubletfinder results column so it is consistent for merging 
for (i in range){
  #pull the name of the df results meta data column 
  df_label <- df_col[[i]]
  #Paste a string to call the metadata column
  meta_col <- paste0("Sce[[i]]@meta.data$",df_label)
  #And paste a command to remove this column for later 
  remove_col <- paste0(meta_col," <- NULL")
  #Generate a temporary list of the doublet identities 
  col <- eval(parse(text=meta_col))
  df_results[[i]] <- col
  #Assign these doublet identities to the doublet column 
  Sce[[i]]$doublet <- col
  #Then at the end, remove the DF column to clean up metadata
  #eval(parse(text=remove_col))
}

saveRDS(df_col, file = paste0(outdir,fig_lab,"_DF_parameters.rds"))

#Then load in original object, add DF column,  and save 
Sce <- list()
for (i in range){
  Sce[[i]] <- readRDS(Files[[i]])
  Sce[[i]] <- subset(Sce[[i]], subset = percent.mt < mt_cutoff)
  Sce[[i]] <- subset(Sce[[i]], subset = nFeature_RNA > nfeature_cutoff)
  Sce[[i]] <- subset(Sce[[i]], subset = nCount_RNA > ncount_cutoff)
  Sce[[i]]$DoubletFinder <- df_results[[i]]
  DefaultAssay(Sce[[i]]) <- "decontXcounts"
  saveRDS(Sce[[i]], file = paste0(outdir,fig_lab,i,"_DF_scater_object.rds"))
}
