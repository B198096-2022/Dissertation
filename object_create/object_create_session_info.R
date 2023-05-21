> sessionInfo()
R version 4.2.3 (2023-03-15)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.6

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] glmGamPoi_1.9.0             ggplot2_3.4.2              
 [3] DoubletFinder_2.0.3         sctransform_0.3.5          
 [5] patchwork_1.1.2             dplyr_1.1.1                
 [7] Seurat_4.3.0                SeuratObject_4.1.3         
 [9] sp_1.6-0                    celda_1.14.2               
[11] Matrix_1.5-4                SingleCellExperiment_1.20.1
[13] SummarizedExperiment_1.28.0 Biobase_2.58.0             
[15] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
[17] IRanges_2.32.0              S4Vectors_0.36.2           
[19] BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
[21] matrixStats_0.63.0         

loaded via a namespace (and not attached):
  [1] utf8_1.2.3                 spatstat.explore_3.1-0     reticulate_1.28           
  [4] R.utils_2.12.2             tidyselect_1.2.0           htmlwidgets_1.6.2         
  [7] grid_4.2.3                 combinat_0.0-8             BiocParallel_1.32.6       
 [10] Rtsne_0.16                 munsell_0.5.0              ScaledMatrix_1.6.0        
 [13] ragg_1.2.5                 codetools_0.2-19           ica_1.0-3                 
 [16] future_1.32.0              miniUI_0.1.1.1             withr_2.5.0               
 [19] spatstat.random_3.1-4      colorspace_2.1-0           progressr_0.13.0          
 [22] knitr_1.42                 rstudioapi_0.14            ROCR_1.0-11               
 [25] assertive.base_0.0-9       tensor_1.5                 listenv_0.9.0             
 [28] labeling_0.4.2             GenomeInfoDbData_1.2.9     polyclip_1.10-4           
 [31] farver_2.1.1               parallelly_1.35.0          vctrs_0.6.2               
 [34] generics_0.1.3             xfun_0.38                  R6_2.5.1                  
 [37] doParallel_1.0.17          ggbeeswarm_0.7.1           rsvd_1.0.5                
 [40] RcppEigen_0.3.3.9.3        bitops_1.0-7               spatstat.utils_3.0-2      
 [43] gridGraphics_0.5-1         DelayedArray_0.24.0        promises_1.2.0.1          
 [46] scales_1.2.1               beeswarm_0.4.0             gtable_0.3.3              
 [49] beachmat_2.14.2            globals_0.16.2             goftest_1.2-3             
 [52] spam_2.9-1                 rlang_1.1.0                systemfonts_1.0.4         
 [55] splines_4.2.3              lazyeval_0.2.2             spatstat.geom_3.1-0       
 [58] yaml_2.3.7                 reshape2_1.4.4             abind_1.4-5               
 [61] httpuv_1.6.9               tools_4.2.3                ellipsis_0.3.2            
 [64] RColorBrewer_1.1-3         ggridges_0.5.4             Rcpp_1.0.10               
 [67] plyr_1.8.8                 sparseMatrixStats_1.10.0   zlibbioc_1.44.0           
 [70] purrr_1.0.1                RCurl_1.98-1.12            dbscan_1.1-11             
 [73] deldir_1.0-6               pbapply_1.7-0              viridis_0.6.2             
 [76] cowplot_1.1.1              zoo_1.8-12                 ggrepel_0.9.3             
 [79] cluster_2.1.4              magrittr_2.0.3             data.table_1.14.8         
 [82] magick_2.7.4               scattermore_0.8            lmtest_0.9-40             
 [85] RANN_2.6.1                 fitdistrplus_1.1-11        mime_0.12                 
 [88] evaluate_0.20              xtable_1.8-4               gridExtra_2.3             
 [91] compiler_4.2.3             scater_1.26.1              tibble_3.2.1              
 [94] KernSmooth_2.23-20         R.oo_1.25.0                htmltools_0.5.5           
 [97] later_1.3.1                tidyr_1.3.0                MCMCprecision_0.4.0       
[100] assertive.files_0.0-2      MASS_7.3-58.3              assertive.numbers_0.0-2   
[103] cli_3.6.1                  assertive.types_0.0-3      R.methodsS3_1.8.2         
[106] parallel_4.2.3             dotCall64_1.0-2            igraph_1.4.2              
[109] pkgconfig_2.0.3            plotly_4.10.1              scuttle_1.8.4             
[112] spatstat.sparse_3.0-1      foreach_1.5.2              vipor_0.4.5               
[115] XVector_0.38.0             stringr_1.5.0              digest_0.6.31             
[118] RcppAnnoy_0.0.20           spatstat.data_3.0-1        rmarkdown_2.21            
[121] leiden_0.4.3               enrichR_3.1                uwot_0.1.14               
[124] DelayedMatrixStats_1.20.0  curl_5.0.0                 shiny_1.7.4               
[127] rjson_0.2.21               lifecycle_1.0.3            nlme_3.1-162              
[130] jsonlite_1.8.4             BiocNeighbors_1.16.0       viridisLite_0.4.2         
[133] fansi_1.0.4                pillar_1.9.0               lattice_0.21-8            
[136] ggrastr_1.0.1              fastmap_1.1.1              httr_1.4.5                
[139] survival_3.5-5             glue_1.6.2                 FNN_1.1.3.2               
[142] png_0.1-8                  iterators_1.0.14           multipanelfigure_2.1.2    
[145] assertive.properties_0.0-5 stringi_1.7.12             textshaping_0.3.6         
[148] BiocSingular_1.14.0        irlba_2.3.5.1              future.apply_1.10.0 
