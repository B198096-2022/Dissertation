
#setup libraries and cwd
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

#BiocManager::install("glmGamPoi")
#install.packages("glmGamPoi")
#library(glmGamPoi)


setwd("/Volumes/cmvm/scs/groups/pramacha-GROUP/Max_Hammer")
######################################################################################
sce <- Read10X("./cellranger/zhang/zhang_healthy_1/outs/filtered_feature_bc_matrix/")
zh1 <- CreateSeuratObject(counts = sce, project = "zh1", min.cells = 3, min.features = 200)
zh1

#add column for percent mitochondrial 
zh1[["percent.mt"]] <- PercentageFeatureSet(zh1, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(zh1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

setwd("/Users/maxhammer/Desktop/Diss_QC")

png("zh1_QC")
VlnPlot(zh1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

plot1 <- FeatureScatter(zh1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zh1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

png("zh1_QC2")
plot1 + plot2
dev.off()

#Need to re-upload the data if I want to re-filter (or save it as filtered object)

#Filter the data
#Set nFeatures to > 800 and nCount > 300 to remove debris, percent.mt < 20 to remove damaged cells  
#And cap of nFeatures < 4000 and nCount < 15000 to remove concerningh igh values
zh1_f1 <- subset(zh1, subset = nFeature_RNA > 800 & nFeature_RNA < 4000 & nCount_RNA > 300 & nCount_RNA < 15000 & percent.mt < 15)
plot1 <- FeatureScatter(zh1_f1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zh1_f1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

png("zh1_filtered_QC2")
plot1 + plot2
dev.off()

png("zh1_filtered_QC")
VlnPlot(zh1_f1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

setwd("/Volumes/cmvm/scs/groups/pramacha-GROUP/Max_Hammer")
#Normalize the data 

#I can't get glmGamPoi to work, it says not compatible with this version 
#But the website says compatible with 4.2 and this is 4.2.3 so idk
#zh1_n <- SCTransform(zh1, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

zh1_n <- SCTransform(zh1_f1, vars.to.regress = "percent.mt", verbose = FALSE)

#Run PCA
zh1_p <- RunPCA(zh1_n, verbose = FALSE)

#Run UMAP 
zh1_u <- RunUMAP(zh1_p, dims = 1:30, verbose = FALSE)

#Clustering 
zh1_uc <- FindNeighbors(zh1_u, dims = 1:30, verbose = FALSE)
zh1_uc <- FindClusters(zh1_uc, verbose = FALSE)

DimPlot(zh1_uc, label = TRUE) + NoLegend()

setwd("/Users/maxhammer/Desktop/Diss_QC")

png("zh1_UMAP")
DimPlot(zh1_uc, label = TRUE) + NoLegend()
dev.off()

setwd("/Volumes/cmvm/scs/groups/pramacha-GROUP/Max_Hammer")

#Save the UMAP output to the DataStore

saveRDS(zh1_uc, file = "./analysis/zhang_healthy_1_UMAP.rds")


################################
#    Now for another sample!   #
################################

sce <- Read10X("./cellranger/zhang/zhang_healthy_2/outs/filtered_feature_bc_matrix/")
zh2 <- CreateSeuratObject(counts = sce, project = "zh2", min.cells = 3, min.features = 200)


#add column for percent mitochondrial 
zh2[["percent.mt"]] <- PercentageFeatureSet(zh2, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(zh2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

setwd("/Users/maxhammer/Desktop/Diss_QC")

png("zh2_QC")
VlnPlot(zh2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

plot1 <- FeatureScatter(zh2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zh2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

png("zh2_QC2")
plot1 + plot2
dev.off()

#Filter the data
#Set nFeatures to > 800 and nCount > 300 to remove debris, percent.mt < 20 to remove damaged cells  
#And cap of nFeatures < 4000 and nCount < 15000 to remove concerningh igh values
zh2_f1 <- subset(zh2, subset = nFeature_RNA > 800 & nFeature_RNA < 4000 & nCount_RNA > 300 & nCount_RNA < 15000 & percent.mt < 15)
plot1 <- FeatureScatter(zh2_f1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zh2_f1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

png("zh2_filtered_QC2")
plot1 + plot2
dev.off()

png("zh2_filtered_QC")
VlnPlot(zh2_f1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

setwd("/Volumes/cmvm/scs/groups/pramacha-GROUP/Max_Hammer")

#Normalize the data 


zh2_n <- SCTransform(zh2_f1, vars.to.regress = "percent.mt", verbose = FALSE)

#Run PCA
zh2_p <- RunPCA(zh2_n, verbose = FALSE)

#Run UMAP 
zh2_u <- RunUMAP(zh2_p, dims = 1:30, verbose = FALSE)

#Clustering 
zh2_uc <- FindNeighbors(zh2_u, dims = 1:30, verbose = FALSE)
zh2_uc <- FindClusters(zh2_uc, verbose = FALSE)

DimPlot(zh2_uc, label = TRUE) + NoLegend()

setwd("/Users/maxhammer/Desktop/Diss_QC")

png("zh2_UMAP")
DimPlot(zh2_uc, label = TRUE) + NoLegend()
dev.off()

setwd("/Volumes/cmvm/scs/groups/pramacha-GROUP/Max_Hammer")

#Save the UMPA output to the DataStore

saveRDS(zh2_uc, file = "./analysis/zhang_healthy_2_UMAP.rds")


################################
#    And another sample!       #
################################

sce <- Read10X("./cellranger/zhang/zhang_healthy_3/outs/filtered_feature_bc_matrix/")
zh3 <- CreateSeuratObject(counts = sce, project = "zh3", min.cells = 3, min.features = 200)


#add column for percent mitochondrial 
zh3[["percent.mt"]] <- PercentageFeatureSet(zh3, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(zh3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

setwd("/Users/maxhammer/Desktop/Diss_QC")

png("zh3_QC")
VlnPlot(zh3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

plot1 <- FeatureScatter(zh3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zh3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

png("zh3_QC2")
plot1 + plot2
dev.off()


#Filter the data
#Set nFeatures to > 800 and nCount > 300 to remove debris, percent.mt < 20 to remove damaged cells  
#And cap of nFeatures < 4000 and nCount < 15000 to remove concerningh igh values
zh3_f1 <- subset(zh3, subset = nFeature_RNA > 800 & nFeature_RNA < 4000 & nCount_RNA > 300 & nCount_RNA < 15000 & percent.mt < 15)
plot1 <- FeatureScatter(zh3_f1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zh3_f1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

png("zh3_filtered_QC2")
plot1 + plot2
dev.off()

png("zh3_filtered_QC")
VlnPlot(zh3_f1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

setwd("/Volumes/cmvm/scs/groups/pramacha-GROUP/Max_Hammer")


setwd("/Volumes/cmvm/scs/groups/pramacha-GROUP/Max_Hammer")

#Normalize the data 


zh3_n <- SCTransform(zh3_f1, vars.to.regress = "percent.mt", verbose = FALSE)

#Run PCA
zh3_p <- RunPCA(zh3_n, verbose = FALSE)

#Run UMAP 
zh3_u <- RunUMAP(zh3_p, dims = 1:30, verbose = FALSE)

#Clustering 
zh3_uc <- FindNeighbors(zh3_u, dims = 1:30, verbose = FALSE)
zh3_uc <- FindClusters(zh3_uc, verbose = FALSE)

DimPlot(zh3_uc, label = TRUE) + NoLegend()

#Save the UMPA output to the DataStore




