Source codes for Fig.1

sessionInfo()
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: CentOS Linux 8

Matrix products: default
BLAS:   /data3/software/R/4.4.1/lib64/R/lib/libRblas.so 
LAPACK: /data3/software/R/4.4.1/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Asia/Shanghai
tzcode source: system (glibc)

attached base packages:
 [1] compiler  grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RColorBrewer_1.1-3          reshape2_1.4.4              ggthemes_5.1.0             
 [4] tidydr_0.0.5                ggsci_3.2.0                 circlize_0.4.16            
 [7] trend_1.1.6                 gtools_3.9.5                clusterProfiler_4.14.6     
[10] ComplexHeatmap_2.22.0       cowplot_1.1.3               ImpulseDE2_0.99.10         
[13] ggprism_1.0.5               ggpubr_0.6.0                gghalves_0.1.4             
[16] msigdbr_10.0.2              AUCell_1.28.0               scCustomize_3.0.1          
[19] data.table_1.17.2           mascarade_0.2.1             SeuratData_0.2.2.9002      
[22] Augur_1.0.3                 viridis_0.6.5               viridisLite_0.4.2          
[25] ggbeeswarm_0.7.2            readxl_1.4.5                lubridate_1.9.4            
[28] forcats_1.0.0               purrr_1.0.4                 readr_2.1.5                
[31] tidyr_1.3.1                 tibble_3.2.1                tidyverse_2.0.0            
[34] BiocParallel_1.40.2         qs_0.27.3                   patchwork_1.3.0            
[37] dplyr_1.1.4                 scran_1.34.0                scater_1.34.1              
[40] ggplot2_3.5.2               scuttle_1.16.0              SingleCellExperiment_1.28.1
[43] SummarizedExperiment_1.36.0 Biobase_2.66.0              GenomicRanges_1.58.0       
[46] GenomeInfoDb_1.42.3         IRanges_2.40.1              S4Vectors_0.44.0           
[49] BiocGenerics_0.52.0         MatrixGenerics_1.21.0       matrixStats_1.5.0          
[52] miloR_2.2.0                 edgeR_4.4.2                 limma_3.62.2               
[55] Seurat_5.3.0                SeuratObject_5.1.0          sp_2.2-0                   
[58] stringr_1.5.1              

loaded via a namespace (and not attached):
  [1] igraph_2.1.4              graph_1.84.1              ica_1.0-3                 plotly_4.10.4            
  [5] Formula_1.2-5             rematch2_2.1.2            maps_3.4.2.1              devtools_2.4.5           
  [9] zlibbioc_1.52.0           tidyselect_1.2.1          bit_4.6.0                 doParallel_1.0.17        
 [13] clue_0.3-66               lattice_0.22-7            rjson_0.2.23              blob_1.2.4               
 [17] urlchecker_1.0.1          S4Arrays_1.6.0            parallel_4.4.1            dichromat_2.0-0.1        
 [21] png_0.1-8                 ggplotify_0.1.2           cli_3.6.5                 goftest_1.2-3            
 [25] textshaping_1.0.1         pbmcapply_1.5.1           bluster_1.12.0            BiocNeighbors_2.0.1      
 [29] uwot_0.2.3                curl_6.2.2                tidytree_0.4.6            mime_0.13                
 [33] evaluate_1.0.3            stringi_1.8.7             backports_1.5.0           desc_1.4.3               
 [37] XML_3.99-0.18             httpuv_1.6.16             AnnotationDbi_1.68.0      paletteer_1.6.0          
 [41] magrittr_2.0.3            rappdirs_0.3.3            splines_4.4.1             prodlim_2025.04.28       
 [45] RApiSerialize_0.1.4       ggraph_2.2.1              sctransform_0.4.2         sessioninfo_1.2.3        
 [49] DBI_1.2.3                 withr_3.0.2               class_7.3-23              systemfonts_1.2.3        
 [53] enrichplot_1.26.6         lmtest_0.9-40             GSEABase_1.68.0           tidygraph_1.3.1          
 [57] BiocManager_1.30.26       htmlwidgets_1.6.4         fs_1.6.6                  biomaRt_2.62.1           
 [61] ggrepel_0.9.6             labeling_0.4.3            SparseArray_1.6.2         DESeq2_1.48.1            
 [65] cellranger_1.1.0          annotate_1.84.0           reticulate_1.42.0         zoo_1.8-14               
 [69] XVector_0.46.0            knitr_1.50                UCSC.utils_1.2.0          timechange_0.3.0         
 [73] foreach_1.5.2             ggtree_3.14.0             timeDate_4041.110         R.oo_1.27.1              
 [77] RSpectra_0.16-2           irlba_2.3.5.1             tester_0.2.0              ggrastr_1.0.2            
 [81] gridGraphics_0.5-1        fastDummies_1.7.5         ellipsis_0.3.2            lazyeval_0.2.2           
 [85] survival_3.8-3            scattermore_1.2           crayon_1.5.3              RcppAnnoy_0.0.22         
 [89] progressr_0.15.1          tweenr_2.0.3              mapproj_1.2.11            later_1.4.2              
 [93] ggridges_0.5.6            codetools_0.2-20          GlobalOptions_0.1.2       profvis_0.4.0            
 [97] KEGGREST_1.46.0           Rtsne_0.17                shape_1.4.6.1             filelock_1.0.3           
[101] pkgconfig_2.0.3           xml2_1.3.8                spatstat.univar_3.1-3     aplot_0.2.5              
[105] ape_5.8-1                 spatstat.sparse_3.1-0     xtable_1.8-4              extraDistr_1.10.0        
[109] car_3.1-3                 plyr_1.8.9                httr_1.4.7                tools_4.4.1              
[113] globals_0.18.0            hardhat_1.4.1             pkgbuild_1.4.7            beeswarm_0.4.0           
[117] broom_1.0.8               nlme_3.1-168              dbplyr_2.5.0              assertthat_0.2.1         
[121] digest_0.6.37             numDeriv_2016.8-1.1       Matrix_1.7-3              furrr_0.3.1              
[125] farver_2.1.2              tzdb_0.5.0                yulab.utils_0.2.0         rpart_4.1.24             
[129] glue_1.8.0                cachem_1.1.0              BiocFileCache_2.14.0      polyclip_1.10-7          
[133] generics_0.1.4            Biostrings_2.74.1         rsample_1.3.0             parallelly_1.44.0        
[137] pkgload_1.4.0             statmod_1.5.0             RcppHNSW_0.6.0            ragg_1.4.0               
[141] ScaledMatrix_1.14.0       carData_3.0-5             pbapply_1.7-2             httr2_1.1.2              
[145] spam_2.11-1               gson_0.1.0                dqrng_0.4.1               utf8_1.2.5               
[149] gower_1.0.2               graphlayouts_1.2.2        ggsignif_0.6.4            gridExtra_2.3            
[153] shiny_1.10.0              lava_1.8.1                GenomeInfoDbData_1.2.13   R.utils_2.13.0           
[157] pals_1.10                 memoise_2.0.1             scales_1.4.0              R.methodsS3_1.8.2        
[161] future_1.49.0             RANN_2.6.2                stringfish_0.16.0         spatstat.data_3.1-6      
[165] rstudioapi_0.17.1         cluster_2.1.8.1           janitor_2.2.1             spatstat.utils_3.1-3     
[169] hms_1.1.3                 fitdistrplus_1.2-2        colorspace_2.1-1          rlang_1.1.6              
[173] DelayedMatrixStats_1.28.1 sparseMatrixStats_1.19.0  ipred_0.9-15              dotCall64_1.2            
[177] ggtangle_0.0.6            ggforce_0.4.2             xfun_0.52                 remotes_2.5.0            
[181] recipes_1.3.0             iterators_1.0.14          abind_1.4-8               randomForest_4.7-1.2     
[185] GOSemSim_2.32.0           treeio_1.30.0             ps_1.9.1                  promises_1.3.2           
[189] qvalue_2.38.0             RSQLite_2.3.11            fgsea_1.32.4              DelayedArray_0.32.0      
[193] GO.db_3.20.0              prettyunits_1.2.0         beachmat_2.22.0           listenv_0.9.1            
[197] Rcpp_1.0.14               parsnip_1.3.2             BiocSingular_1.22.0       tensor_1.5               
[201] usethis_3.1.0             MASS_7.3-65               progress_1.2.3            babelgene_22.9           
[205] spatstat.random_3.3-3     R6_2.6.1                  fastmatch_1.1-6           fastmap_1.2.0            
[209] rstatix_0.7.2             vipor_0.4.7               ROCR_1.0-11               rsvd_1.0.5               
[213] nnet_7.3-20               gtable_0.3.6              KernSmooth_2.23-26        miniUI_0.1.2             
[217] deldir_2.0-4              htmltools_0.5.8.1         yardstick_1.3.2           RcppParallel_5.1.10      
[221] bit64_4.6.0-1             spatstat.explore_3.4-2    lifecycle_1.0.4           msigdbdf_24.1.1          
[225] processx_3.8.6            callr_3.7.6               vctrs_0.6.5               spatstat.geom_3.3-6      
[229] snakecase_0.11.1          DOSE_4.0.1                ggfun_0.1.8               future.apply_1.11.3      
[233] pracma_2.4.4              pillar_1.10.2             metapod_1.14.0            locfit_1.5-9.12          
[237] jsonlite_2.0.0            GetoptLong_1.0.5 
