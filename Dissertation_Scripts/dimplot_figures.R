library(Seurat)
library(ggplot2)


all_merged <- readRDS("~/Desktop/big_pilot/merged_objects/pilot2_object_scVI_harmony.rds")
small_merged <- readRDS("~/Desktop/big_pilot/merged_objects/pilot2_object_rpca_cca.rds")


tamburini_all <- subset(all_merged, subset = (author == "tamburini"))
guilliams_all <- subset(all_merged, subset = (author == "guilliams"))
zhang_all <- subset(all_merged, subset = (author == "zhang"))
ramachandran_all <- subset(all_merged, subset = (author == "ramachandran"))
buonomo_all <- subset(all_merged, subset = (author == "buonomo"))

tamburini_small <- subset(small_merged, subset = (author == "tamburini"))
guilliams_small <- subset(small_merged, subset = (author == "guilliams"))
zhang_small <- subset(small_merged, subset = (author == "zhang"))
ramachandran_small <- subset(small_merged, subset = (author == "ramachandran"))
buonomo_small <- subset(small_merged, subset = (author == "buonomo"))

fibrosis_n <- subset(all_merged, subset = (fibrosis_stage == "none-low"))
fibrosis_m <- subset(all_merged, subset = (fibrosis_stage == "moderate"))
fibrosis_s <- subset(all_merged, subset = (fibrosis_stage == "severe"))




int_plt_harmony_patient <- DimPlot(
  all_merged,
  reduction = "umap.harmony",
  group.by = c("patient"))

int_plt_scvi_patient <- DimPlot(
  all_merged,
  reduction = "umap.scvi",
  group.by = c("patient"))

int_plt_cca_patient <- DimPlot(
  small_merged,
  reduction = "umap.cca",
  group.by = c("patient"))

int_plt_rpca_patient <- DimPlot(
  small_merged,
  reduction = "umap.rpca",
  group.by = c("patient"))
int_plt_rpca_patient + theme(text = element_text(size = 25))


unint_plt_patient <- DimPlot(
  all_merged,
  reduction = "umap.unintegrated",
  group.by = c("patient"))

unint_plt_patient + theme(text = element_text(size = 25))
int_plt_harmony_patient + theme(text = element_text(size = 25))
int_plt_cca_patient + theme(text = element_text(size = 25))
int_plt_rpca_patient + theme(text = element_text(size = 25))
int_plt_scvi_patient + theme(text = element_text(size = 25))




harmony_cell_types <- DimPlot(
  all_merged,
  reduction = "umap.harmony",
  group.by = c("cell_type"), 
  label = TRUE,
  cols = c("#f8766d", "#e58600", "#a3a500", "#cccccc", "#6e6d6d" ,"#00ba39", "#00bf7d", "#00c0af", "#01bcd7", 
           "#00b0f6", "#619cff", "#b983ff", "#e76bf3", "#fd60d1", "#ff68a4")
)

unint_cell_types <- DimPlot(
  all_merged,
  reduction = "umap.unintegrated",
  group.by = c("cell_type"), 
  label = TRUE,
  cols = c("#f8766d", "#e58600", "#a3a500", "#cccccc", "#6e6d6d" ,"#00ba39", "#00bf7d", "#00c0af", "#01bcd7", 
           "#00b0f6", "#619cff", "#b983ff", "#e76bf3", "#fd60d1", "#ff68a4")
)

scvi_cell_types <- DimPlot(
  all_merged,
  reduction = "umap.scvi",
  group.by = c("cell_type"), 
  label = TRUE,
  cols = c("#f8766d", "#e58600", "#a3a500", "#cccccc", "#6e6d6d" ,"#00ba39", "#00bf7d", "#00c0af", "#01bcd7", 
           "#00b0f6", "#619cff", "#b983ff", "#e76bf3", "#fd60d1", "#ff68a4")
)

rpca_cell_types <- DimPlot(
  small_merged,
  reduction = "umap.rpca",
  group.by = c("cell_type"), 
  label = TRUE,
  cols = c("#f8766d", "#e58600", "#a3a500", "#cccccc", "#6e6d6d" ,"#00ba39", "#00bf7d", "#00c0af", "#01bcd7", 
           "#00b0f6", "#619cff", "#b983ff", "#e76bf3", "#fd60d1", "#ff68a4")
)

cca_cell_types <- DimPlot(
  small_merged,
  reduction = "umap.cca",
  group.by = c("cell_type"), 
  label = TRUE,
  cols = c("#f8766d", "#e58600", "#a3a500", "#cccccc", "#6e6d6d" ,"#00ba39", "#00bf7d", "#00c0af", "#01bcd7", 
           "#00b0f6", "#619cff", "#b983ff", "#e76bf3", "#fd60d1", "#ff68a4")
)

harmony_cell_types + theme(text = element_text(size = 20))
scvi_cell_types + theme(text = element_text(size = 20))
unint_cell_types + theme(text = element_text(size = 20))
rpca_cell_types + theme(text = element_text(size = 20))
cca_cell_types + theme(text = element_text(size = 20))

sm_cell_col <- readRDS("~/Desktop/big_pilot/cell_types/second/small_cell_types_col.rds")
small_merged$cell_type <- sm_cell_col


################
#Fibrosis stage overlay 
int_plt_harmony_fibrosis <- DimPlot(
  all_merged,
  reduction = "umap.harmony",
  group.by = c("fibrosis_stage"),
  cols = c("#e58600", "#00ba39", "#fb4336"))

unint_plt_fibrosis <- DimPlot(
  all_merged,
  reduction = "umap.unintegrated",
  group.by = c("fibrosis_stage"),
  cols = c("#e58600", "#00ba39", "#fb4336"))

unint_plt_fibrosis + theme(text = element_text(size = 25))
int_plt_harmony_fibrosis + theme(text = element_text(size = 25))



################
#Fibrosis stage individual 
int_plt_harmony_fibrosis_n <- DimPlot(
  fibrosis_n,
  reduction = "umap.harmony",
  group.by = c("fibrosis_stage"),
  cols = c("#00ba39"))

int_plt_harmony_fibrosis_n + theme(text = element_text(size = 25))

int_plt_harmony_fibrosis_m <- DimPlot(
  fibrosis_m,
  reduction = "umap.harmony",
  group.by = c("fibrosis_stage"),
  cols = c("#e58600"))

int_plt_harmony_fibrosis_m + theme(text = element_text(size = 25))


int_plt_harmony_fibrosis_s <- DimPlot(
  fibrosis_s,
  reduction = "umap.harmony",
  group.by = c("fibrosis_stage"),
  cols = c("#fb4336"))

int_plt_harmony_fibrosis_s + theme(text = element_text(size = 25))

library(ggpubr)
ggarrange(int_plt_harmony_fibrosis_n, int_plt_harmony_fibrosis_m, int_plt_harmony_fibrosis_s,
          labels = c("None-Low", "Moderate", "Severe"),
          ncol = 1, nrow = 3)



#Unintegrated

unint_plt_fibrosis_n <- DimPlot(
  fibrosis_n,
  reduction = "umap.unintegrated",
  group.by = c("fibrosis_stage"),
  cols = c("#00ba39"))

unint_plt_fibrosis_n + theme(text = element_text(size = 25))

unint_plt_fibrosis_m <- DimPlot(
  fibrosis_m,
  reduction = "umap.unintegrated",
  group.by = c("fibrosis_stage"),
  cols = c("#e58600"))

unint_plt_fibrosis_m + theme(text = element_text(size = 25))


unint_plt_fibrosis_s <- DimPlot(
  fibrosis_s,
  reduction = "umap.unintegrated",
  group.by = c("fibrosis_stage"),
  cols = c("#fb4336"))

unint_plt_fibrosis_s + theme(text = element_text(size = 25))

library(ggpubr)
ggarrange(unint_plt_fibrosis_n, unint_plt_fibrosis_m, unint_plt_fibrosis_s,
          labels = c("None-Low", "Moderate", "Severe"),
          ncol = 1, nrow = 3)

########################
#Study overlay 

int_plt_harmony_study <- DimPlot(
  all_merged,
  reduction = "umap.harmony",
  group.by = c("author"))

unint_plt_study <- DimPlot(
  all_merged,
  reduction = "umap.unintegrated",
  group.by = c("author"))

int_plt_harmony_study + theme(text = element_text(size = 25))
unint_plt_study + theme(text = element_text(size = 25))


####################
#study individual 

unint_plt_study_tamburini <- DimPlot(
  tamburini_all,
  reduction = "umap.unintegrated",
  group.by = c("author"),
  cols = c("#00b0f6"))

unint_plt_study_tamburini + theme(text = element_text(size = 25))

unint_plt_study_buonomo <- DimPlot(
  buonomo_all,
  reduction = "umap.unintegrated",
  group.by = c("author"),
  cols = c("#f8766d"))

unint_plt_study_buonomo + theme(text = element_text(size = 25))

unint_plt_study_guilliams <- DimPlot(
  guilliams_all,
  reduction = "umap.unintegrated",
  group.by = c("author"),
  cols = c("#a3a500"))

unint_plt_study_guilliams + theme(text = element_text(size = 25))


unint_plt_study_ramachandran <- DimPlot(
  ramachandran_all,
  reduction = "umap.unintegrated",
  group.by = c("author"),
  cols = c("#00bf7d"))

unint_plt_study_ramachandran + theme(text = element_text(size = 25))


unint_plt_study_zhang <- DimPlot(
  zhang_all,
  reduction = "umap.unintegrated",
  group.by = c("author"),
  cols = c("#b983ff"))

unint_plt_study_zhang + theme(text = element_text(size = 25))


#Harmony

harmony_plt_study_tamburini <- DimPlot(
  tamburini_all,
  reduction = "umap.harmony",
  group.by = c("author"),
  cols = c("#00b0f6"))

harmony_plt_study_tamburini + theme(text = element_text(size = 25))

harmony_plt_study_buonomo <- DimPlot(
  buonomo_all,
  reduction = "umap.harmony",
  group.by = c("author"),
  cols = c("#f8766d"))

harmony_plt_study_buonomo + theme(text = element_text(size = 25))

harmony_plt_study_guilliams <- DimPlot(
  guilliams_all,
  reduction = "umap.harmony",
  group.by = c("author"),
  cols = c("#a3a500"))

harmony_plt_study_guilliams + theme(text = element_text(size = 25))


harmony_plt_study_ramachandran <- DimPlot(
  ramachandran_all,
  reduction = "umap.harmony",
  group.by = c("author"),
  cols = c("#00bf7d"))

harmony_plt_study_ramachandran + theme(text = element_text(size = 25))


harmony_plt_study_zhang <- DimPlot(
  zhang_all,
  reduction = "umap.harmony",
  group.by = c("author"),
  cols = c("#b983ff"))

harmony_plt_study_zhang + theme(text = element_text(size = 25))


###### QC Plots
#cols = c("#05bfc4", "#f7766d")

harmony_db <- DimPlot(
  all_merged,
  reduction = "umap.harmony",
  group.by = c("DoubletFinder"),
  cols = c("#f7766d", "#cccccc"))

harmony_db + theme(text = element_text(size = 25))

harmony_mo <- DimPlot(
  all_merged,
  reduction = "umap.harmony",
  group.by = c("mt_outlier"),
  cols = c("#cccccc", "#f7766d"))

harmony_mo + theme(text = element_text(size = 25))

harmony_co <- DimPlot(
  all_merged,
  reduction = "umap.harmony",
  group.by = c("count_outlier"),
  cols = c("#cccccc", "#f7766d"))

harmony_co + theme(text = element_text(size = 25))

harmony_fo <- DimPlot(
  all_merged,
  reduction = "umap.harmony",
  group.by = c("detected_outlier"),
  cols = c("#cccccc", "#f7766d"))

harmony_fo + theme(text = element_text(size = 25))

harmony_to <- DimPlot(
  all_merged,
  reduction = "umap.harmony",
  group.by = c("final_outlier"), 
  cols = c("#cccccc", "#f7766d"))

harmony_to + theme(text = element_text(size = 25))

#

harmony_clusters <- DimPlot(
  all_merged,
  reduction = "umap.harmony",
  group.by = c("harmony_clusters"),
  label = TRUE)

harmony_clusters + theme(text = element_text(size = 20))

cca_clusters <- DimPlot(
  small_merged,
  reduction = "umap.cca",
  group.by = c("cca_clusters"),
  label = TRUE)

cca_clusters + theme(text = element_text(size = 20))

rpca_clusters <- DimPlot(
  small_merged,
  reduction = "umap.rpca",
  group.by = c("rpca_clusters"),
  label = TRUE)

rpca_clusters + theme(text = element_text(size = 20))

unint_clusters <- DimPlot(
  all_merged,
  reduction = "umap.unintegrated",
  group.by = c("unintegrated_clusters"),
  label = TRUE)

unint_clusters + theme(text = element_text(size = 20))

scvi_clusters <- DimPlot(
  all_merged,
  reduction = "umap.scvi",
  group.by = c("scvi_clusters"),
  label = TRUE)

scvi_clusters + theme(text = element_text(size = 20))


###############################################
#no labels


harmony_clusters <- DimPlot(
  all_merged,
  reduction = "umap.harmony",
  group.by = c("harmony_clusters"),
  label = FALSE)


cca_clusters <- DimPlot(
  small_merged,
  reduction = "umap.cca",
  group.by = c("cca_clusters"),
  label = FALSE)

rpca_clusters <- DimPlot(
  small_merged,
  reduction = "umap.rpca",
  group.by = c("rpca_clusters"),
  label = FALSE)

unint_clusters <- DimPlot(
  all_merged,
  reduction = "umap.unintegrated",
  group.by = c("unintegrated_clusters"),
  label = FALSE)


scvi_clusters <- DimPlot(
  all_merged,
  reduction = "umap.scvi",
  group.by = c("scvi_clusters"),
  label = FALSE)










#
