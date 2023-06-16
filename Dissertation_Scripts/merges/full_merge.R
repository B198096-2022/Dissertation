old_lib <- "/Library/Frameworks/R.framework/Versions/4.2/Resources/old-versions"
library("SeuratObject", lib.loc = old_lib)
library("Seurat", lib.loc = old_lib)


options(future.globals.maxSize = 1e12)
library(dplyr)

datadir <- "~/Desktop/big_pilot/merged_objects/"
#Ramachandran
rhc <- readRDS(paste0(datadir,"rhc_merged_object.rds"))
rfc <- readRDS(paste0(datadir,"rfc_merged_object.rds"))

#Buonomo
bhc <- readRDS(paste0(datadir,"bhc_merged_object.rds"))
bfc <- readRDS(paste0(datadir,"bfc_merged_object.rds"))

#Guilliams
ghc <- readRDS(paste0(datadir,"ghc_merged_object.rds"))
gfc <- readRDS(paste0(datadir,"gfc_merged_object.rds"))

#Tamburini
thc <- readRDS(paste0(datadir,"thc_merged_object.rds"))
tfc <- readRDS(paste0(datadir,"tfc_merged_object.rds"))

#Zhang
zhc <- readRDS(paste0(datadir,"zhc_merged_object.rds"))
zfc <- readRDS(paste0(datadir,"zfc_merged_object.rds"))

full_merged <- merge(rhc, y = c(rfc, bhc, bfc, ghc, gfc, thc, tfc, zhc, zfc), project ="full_merged")
saveRDS(full_merged, file = paste0(datadir, "full_merged_object.rds"))





