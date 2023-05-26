#########################
# ghc4, which is ghc_p3 #
# is labeled gfc_p3     #
#########################

##################################
# The merge includes this error  #
# and also repeats bfc1 twice    #
##################################

#Fixing all of the ghc3 containing objects 
fix1 <- readRDS("~/Desktop/Bioinformatics/Dissertation/guilliams_fix/ghc4_DF_scater_object.rds")
fix1$patient <- "ghc_p3"
saveRDS(fix1, file = "~/Desktop/big_pilot/guilliams/ghc4_DF_scater_object.rds")

fix2 <- readRDS("~/Desktop/Bioinformatics/Dissertation/guilliams_fix/ghc4_scater_object.rds")
fix2$patient <- "ghc_p3"
saveRDS(fix2, file = "~/Desktop/Bioinformatics/Dissertation/guilliams_fix/ghc4_scater_object.rds")

fix3 <- readRDS("~/Desktop/Bioinformatics/Dissertation/guilliams_fix/ghc4_dx_SeuratObject.rds")
fix3$patient <- "ghc_p3"
saveRDS(fix3, file = "~/Desktop/Bioinformatics/Dissertation/guilliams_fix/ghc4_dx_SeuratObject.rds")

fix4 <- readRDS("~/Desktop/Bioinformatics/Dissertation/guilliams_fix/ghc_merged_object.rds")

dim(fix4)

for ( i in 1:dim(fix4)[2]){
  if (as.character(fix4$patient[i]) == "gfc_p3"){
    fix4$patient[i] <- "ghc_p3"
  }
}

saveRDS(fix4, file = "~/Desktop/Bioinformatics/Dissertation/guilliams_fix/ghc_merged_object.rds")

