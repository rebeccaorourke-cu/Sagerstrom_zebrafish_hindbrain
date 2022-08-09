HB13hpf\_DEAB\_neural\_subset
================
Rebecca O’Rourke
8/9/2022

# 1. libraries, palette and functions

``` r
suppressPackageStartupMessages({
  library(Seurat)
  library(org.Dr.eg.db)
  library(BSgenome.Drerio.UCSC.danRer11)
  library(Signac)
  library(knitr)
  library(kableExtra)
  library(dplyr)
  library(ggplot2)
  library(ggsci)
  library(limma)
  library(JASPAR2020)
  library(patchwork)
  library(TFBSTools)
  library(motifmatchr)
  library(harmony)
})
```

    ## Warning: package 'AnnotationDbi' was built under R version 4.1.1

    ## Warning: package 'BiocGenerics' was built under R version 4.1.1

    ## Warning: package 'Biobase' was built under R version 4.1.1

    ## Warning: package 'IRanges' was built under R version 4.1.1

    ## Warning: package 'S4Vectors' was built under R version 4.1.3

    ## Warning: package 'BSgenome' was built under R version 4.1.1

    ## Warning: package 'GenomeInfoDb' was built under R version 4.1.1

    ## Warning: package 'GenomicRanges' was built under R version 4.1.1

    ## Warning: package 'Biostrings' was built under R version 4.1.1

    ## Warning: package 'XVector' was built under R version 4.1.1

    ## Warning: package 'rtracklayer' was built under R version 4.1.1

    ## Warning: package 'limma' was built under R version 4.1.3

    ## Warning: package 'TFBSTools' was built under R version 4.1.1

    ## Warning: package 'motifmatchr' was built under R version 4.1.1

``` r
options(future.globals.maxSize = 4000 * 1024^2)
```

``` r
mypal <- pal_igv(palette = "default",alpha = 1)(50)
```

``` r
GetUMAPandClusters <- function(seurat){
  # RNA analysis
  DefaultAssay(seurat) <- "RNA"
  seurat <- SCTransform(seurat, verbose = FALSE, return.only.var.genes = FALSE) 
  seurat <- RunPCA(seurat) 
  ElbowPlot(seurat,ndims = 50)
  seurat <- RunUMAP(seurat, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  
  # ATAC analysis
  # We exclude the first dimension as this is typically correlated with sequencing depth
  DefaultAssay(seurat) <- "peaks"
  
  seurat <- FindTopFeatures(seurat, min.cutoff = 5)
  seurat <- RunTFIDF(seurat)
  seurat <- RunSVD(seurat)
  seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  
  DefaultAssay(seurat) <- "SCT"
  seurat <- FindMultiModalNeighbors(seurat, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  seurat <- RunUMAP(seurat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", assay = "RNA")
  seurat <- FindClusters(seurat, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  seurat <- FindClusters(seurat, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 8.0)
  
  p1 <- DimPlot(seurat, reduction = "umap.rna", label = TRUE, label.size = 5, repel = TRUE, group.by = "wsnn_res.0.8") + ggtitle("RNA") + scale_color_igv()
  p2 <- DimPlot(seurat, reduction = "umap.atac", label = TRUE, label.size = 5, repel = TRUE, group.by = "wsnn_res.0.8") + ggtitle("ATAC") + scale_color_igv()
  p3 <- DimPlot(seurat, reduction = "wnn.umap", label = TRUE, label.size = 5, repel = TRUE, group.by = "wsnn_res.0.8") + ggtitle("WNN") + scale_color_igv()
  print(p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
  
  return(seurat)
}
```

# 2. Read data

``` r
HB13hpf.DEAB <- readRDS("RDSfiles/HB13hpf_DEAB.clustered.RDS")
DimPlot(HB13hpf.DEAB, reduction = "wnn.umap", label = T, repel = T, group.by = "wsnn_res.3") + scale_color_igv()
```

![](02d_HB13hpf_DEAB_neural_subset_files/figure-gfm/readData-1.png)<!-- -->

# 3. Add cluster ID’s

``` r
Idents(HB13hpf.DEAB) <- "wsnn_res.3"
HB13hpf.DEAB <- RenameIdents(HB13hpf.DEAB,
                             "0" = "neural crest",
                             "1" = "r1/r2/r4",
                             "2" = "MHB.1",
                             "3" = "MHB.2",
                             "4" = "HB-like",
                             "5" = "MHB.3",
                             "6" = "Placode.1",
                             "7" = "Endoderm.1",
                             "8" = "Eye.1",
                             "9" = "SC?",
                             "10" = "Placode.2",
                             "11" = "Kidney?",
                             "12" = "FB/eye",
                             "13" = "Muscle.1",
                             "14" = "r3",
                             "15" = "cilia/SC?",
                             "16" = "Placode.3",
                             "17" = "MHB.1",
                             "18" = "CHB/SC",
                             "19" = "somite ",
                             "20" = "HiRibosomal",
                             "21" = "Endoderm.2",
                             "22" = "Muscle.2",
                             "23" = "Eye.2",
                             "24" = "heart",
                             "25" = "Notochord",
                             "26" = "Ciliated",
                             "27" = "neurons",
                             "28" = "Skeletal",
                             "29" = "Blood",
                             "30" = "Kidney?",
                             "31" = "skin (epidermis)?")
HB13hpf.DEAB$WholeClusters <- Idents(HB13hpf.DEAB)
DimPlot(HB13hpf.DEAB, reduction = "wnn.umap", label = T, repel = T) + scale_color_igv()
```

![](02d_HB13hpf_DEAB_neural_subset_files/figure-gfm/clusterIDs-1.png)<!-- -->

# 4. Subset neural clusters

the following are neural clusters:

1, 2, 3, 4, 5, 8, 9, 12, 14, 15, 17, 18, 23, 26, 27

``` r
Idents(HB13hpf.DEAB) <- "wsnn_res.3"
neural <- subset(HB13hpf.DEAB, idents = c("1","2","3","4","5","8","9","12","14","15","17","18","23","26","27"))
DimPlot(neural, label = T, repel = T, reduction = "wnn.umap") + scale_color_igv()
```

![](02d_HB13hpf_DEAB_neural_subset_files/figure-gfm/subset-1.png)<!-- -->

# 5. Reprocess only neural cells

Starting with only neural cells, rerun SCTransform, PCA, UMAP, SVD,
Multimodal Neighbors and Clustering

``` r
neural <- GetUMAPandClusters(neural)
```

    ## PC_ 1 
    ## Positive:  rnf220a, cdh11, ntn1b, dnah5, cfap57, cfap45, dnah2, armc4, syne1b, ntn1a 
    ##     pcdh1b, mllt3, dnah6, foxj1a, cacna1db, col11a1a, ptch2, cfap58, tctex1d1, dnah7 
    ##     bicc2, zbbx, nfasca, hoxc3a, chl1a, ccdc13, adamts3, robo3, cfap52, pacrg 
    ## Negative:  pax3a, zbtb16a, pax2a, pax7a, pax5, ror1, meis1b, en2b, efna2a, tenm3 
    ##     ndst2b, pax7b, jag2b, zfhx4, kirrel3l, scube2, wnt8b, wnt3, znrf3, gli2a 
    ##     boc, fam189a1, efna5a, zbtb16b, zeb1b, si:ch211-286o17.1, cachd1, tenm4, nr2f5, gli2b 
    ## PC_ 2 
    ## Positive:  lmx1bb, fras1, si:dkey-91m11.5, tenm4, hmcn2.1, zfhx4, tcf7l2, chl1a, fgfrl1a, plxna1a 
    ##     pax2a, traf4a, hs6st1a, tenm3, pou2f2a.1, tmem108, scube3, foxp4, lhx2b, cdon 
    ##     lef1, vasnb, met, meis1b, chsy1, vcana, nfia, crb2a, six7, bahcc1b 
    ## Negative:  greb1l, chl1b, cntfr, pax7b, kirrel3l, adgrl2a, her4.2, adamts3, cdh7a, auts2b 
    ##     plekhh1, cyp26b1, ntn1a, ndst3, ptprfb, tead1b, fstl1b, CT990561.1, pax6a, adgrl1a 
    ##     epha4b, rasal2, her4.2.1, lamb2l, FRMD5, cxcl12a, sema3fb, dtnba, vgll3, her4.1 
    ## PC_ 3 
    ## Positive:  pax5, lmx1bb, en2b, pax2a, pax7a, efna2a, pax3a, zfhx4, pax2b, col4a6 
    ##     znrf3, rfx4, ptprnb, rnf220a, ndst2b, wnt3, lef1, plxna1a, cachd1, col18a1a 
    ##     nkd1, her11, gas7a, en2a, zbtb16a, wnt1, jag2b, gria3b, hs3st3b1b, chl1a 
    ## Negative:  fgfrl1a, sfrp5, si:dkey-91m11.5, dclk2a, epha7, traf4a, lhx2b, bmpr1ba, chsy1, nfia 
    ##     tcf7l2, col14a1a, vcana, bahcc1b, six7, hmcn2.1, vasnb, pax6b, mllt3, sox5 
    ##     fras1, epha4a, pard6gb, nr2f5, grasp, lhx9, vax2, pappab, lmo4a, rab11fip4a 
    ## PC_ 4 
    ## Positive:  cdon, pax3a, klf7a, cdh6, msx1b, nfasca, col15a1b, unc5b, ebf2, lmx1bb 
    ##     ZNF423, sema3fb, col18a1a, myo10l3, nbeaa, fam189a1, lmx1ba, dlb, myt1a, srrm4 
    ##     lama5, ctnna2, afap1, ror1, boc, wnt3, neo1a, efnb1, nrp1a, zic2b 
    ## Negative:  ntn1a, cdh7a, hs3st3b1b, ntn1b, rnf220a, epha4b, mllt3, pappab, efna2a, tcf7l2 
    ##     chl1b, adamts3, hs6st1a, meis1b, robo1, fgfrl1a, tenm4, ncam1a, sfrp5, ephb2b 
    ##     qkib, sema5ba, col28a2a, crb1, tet3, tcf7l1b, prickle2b, vcana, chl1a, ndst2b 
    ## PC_ 5 
    ## Positive:  elavl3, ntn1a, ebf2, myt1a, dlb, rnf220a, nfasca, ebf3a, ebf3a.1, srrm4 
    ##     onecut2, efna2a, dla, ctnna2, myo10l3, inavaa, myt1b, gadd45gb.1, spsb4a, her4.2.1 
    ##     cbfa2t3, celf3a, ppp1r14ab, SIPA1, madd, notch1a, dbn1, BX530077.1, neurod1, cdc42bpab 
    ## Negative:  cdon, pax3a, col15a1b, msx1b, pax7b, boc, lmx1bb, epha7, sema3fb, fndc3ba 
    ##     ZNF423, ror1, si:ch211-285f17.1, wnt3, col18a1a, cgnl1, fstl1b, neo1a, col7a1l, BX908388.3 
    ##     unc5b, gli3, ptprfb, plxdc2, epha4a, efnb1, lmx1ba, id1, si:dkey-49n23.1, bmpr1ab

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 12:35:13 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 12:35:13 Read 1306 rows and found 50 numeric columns

    ## 12:35:13 Using Annoy for neighbor search, n_neighbors = 30

    ## 12:35:13 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 12:35:14 Writing NN index file to temp file /var/folders/85/cmygnxr155gdnf4bc3_30jyh0000gp/T//RtmpM1eCsy/file116bb41af2519
    ## 12:35:14 Searching Annoy index using 1 thread, search_k = 3000
    ## 12:35:14 Annoy recall = 100%
    ## 12:35:15 Commencing smooth kNN distance calibration using 1 thread
    ## 12:35:17 Initializing from normalized Laplacian + noise
    ## 12:35:17 Commencing optimization for 500 epochs, with 47458 positive edges
    ## 12:35:19 Optimization finished
    ## Performing TF-IDF normalization

    ## Warning in RunTFIDF.default(object = GetAssayData(object = object, slot =
    ## "counts"), : Some features contain 0 total counts

    ## Running SVD
    ## Scaling cell embeddings
    ## 12:35:43 UMAP embedding parameters a = 0.9922 b = 1.112
    ## 12:35:43 Read 1306 rows and found 49 numeric columns
    ## 12:35:43 Using Annoy for neighbor search, n_neighbors = 30
    ## 12:35:43 Building Annoy index with metric = cosine, n_trees = 50
    ## 0%   10   20   30   40   50   60   70   80   90   100%
    ## [----|----|----|----|----|----|----|----|----|----|
    ## **************************************************|
    ## 12:35:43 Writing NN index file to temp file /var/folders/85/cmygnxr155gdnf4bc3_30jyh0000gp/T//RtmpM1eCsy/file116bb51acadff
    ## 12:35:43 Searching Annoy index using 1 thread, search_k = 3000
    ## 12:35:43 Annoy recall = 100%
    ## 12:35:44 Commencing smooth kNN distance calibration using 1 thread
    ## 12:35:46 Initializing from normalized Laplacian + noise
    ## 12:35:46 Commencing optimization for 500 epochs, with 53248 positive edges
    ## 12:35:49 Optimization finished
    ## Calculating cell-specific modality weights
    ## Finding 20 nearest neighbors for each modality.
    ## Calculating kernel bandwidths
    ## Finding multimodal neighbors
    ## Constructing multimodal KNN graph
    ## Constructing multimodal SNN graph
    ## 12:35:52 UMAP embedding parameters a = 0.9922 b = 1.112
    ## 12:35:53 Commencing smooth kNN distance calibration using 1 thread
    ## 12:35:55 Initializing from normalized Laplacian + noise
    ## 12:35:55 Commencing optimization for 500 epochs, with 35782 positive edges
    ## 12:35:57 Optimization finished

![](02d_HB13hpf_DEAB_neural_subset_files/figure-gfm/reprocess-1.png)<!-- -->

``` r
saveRDS(neural, file = "RDSfiles/HB13hpf_DEAB_neural.RDS")
```

``` r
p1 <- DimPlot(neural, group.by = "wsnn_res.0.8", label = T, reduction = "wnn.umap") + NoLegend() + scale_color_igv()
p2 <- DimPlot(neural, group.by = "wsnn_res.8", label = T, reduction = "wnn.umap") + NoLegend() + scale_color_igv()
p <- p1 + p2 
p
```

![](02d_HB13hpf_DEAB_neural_subset_files/figure-gfm/plot1-1.png)<!-- -->

``` r
p1 <- DimPlot(neural, group.by = "wsnn_res.8", label = T, reduction = "wnn.umap") + NoLegend() + scale_color_igv()
p2 <- DimPlot(neural, group.by = "WholeClusters", label = T, reduction = "wnn.umap") + NoLegend() + scale_color_igv()
p1 + p2
```

![](02d_HB13hpf_DEAB_neural_subset_files/figure-gfm/plot2-1.png)<!-- -->

# 6. sessionInfo

``` r
sessionInfo()
```

    ## R version 4.1.0 (2021-05-18)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] harmony_0.1.0                       Rcpp_1.0.7                         
    ##  [3] motifmatchr_1.16.0                  TFBSTools_1.32.0                   
    ##  [5] patchwork_1.1.1                     JASPAR2020_0.99.10                 
    ##  [7] limma_3.50.3                        ggsci_2.9                          
    ##  [9] ggplot2_3.3.5                       dplyr_1.0.7                        
    ## [11] kableExtra_1.3.4                    knitr_1.36                         
    ## [13] Signac_1.2.1                        BSgenome.Drerio.UCSC.danRer11_1.4.2
    ## [15] BSgenome_1.62.0                     rtracklayer_1.54.0                 
    ## [17] Biostrings_2.62.0                   XVector_0.34.0                     
    ## [19] GenomicRanges_1.46.0                GenomeInfoDb_1.30.0                
    ## [21] org.Dr.eg.db_3.14.0                 AnnotationDbi_1.56.1               
    ## [23] IRanges_2.28.0                      S4Vectors_0.32.4                   
    ## [25] Biobase_2.54.0                      BiocGenerics_0.40.0                
    ## [27] SeuratObject_4.0.4                  Seurat_4.0.1                       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] utf8_1.2.2                  R.utils_2.11.0             
    ##   [3] reticulate_1.22             tidyselect_1.1.1           
    ##   [5] poweRlaw_0.70.6             RSQLite_2.2.8              
    ##   [7] htmlwidgets_1.5.4           grid_4.1.0                 
    ##   [9] docopt_0.7.1                BiocParallel_1.28.0        
    ##  [11] Rtsne_0.15                  munsell_0.5.0              
    ##  [13] codetools_0.2-18            ica_1.0-2                  
    ##  [15] future_1.26.1               miniUI_0.1.1.1             
    ##  [17] withr_2.4.2                 colorspace_2.0-2           
    ##  [19] highr_0.9                   rstudioapi_0.13            
    ##  [21] ROCR_1.0-11                 tensor_1.5                 
    ##  [23] listenv_0.8.0               labeling_0.4.2             
    ##  [25] MatrixGenerics_1.6.0        slam_0.1-48                
    ##  [27] GenomeInfoDbData_1.2.7      polyclip_1.10-0            
    ##  [29] bit64_4.0.5                 farver_2.1.0               
    ##  [31] parallelly_1.32.0           vctrs_0.4.1                
    ##  [33] generics_0.1.1              xfun_0.27                  
    ##  [35] lsa_0.73.2                  ggseqlogo_0.1              
    ##  [37] R6_2.5.1                    bitops_1.0-7               
    ##  [39] spatstat.utils_2.2-0        cachem_1.0.6               
    ##  [41] DelayedArray_0.20.0         assertthat_0.2.1           
    ##  [43] promises_1.2.0.1            BiocIO_1.4.0               
    ##  [45] scales_1.1.1                gtable_0.3.0               
    ##  [47] globals_0.15.1              goftest_1.2-3              
    ##  [49] seqLogo_1.60.0              rlang_1.0.3                
    ##  [51] systemfonts_1.0.4           RcppRoll_0.3.0             
    ##  [53] splines_4.1.0               lazyeval_0.2.2             
    ##  [55] spatstat.geom_2.3-0         yaml_2.2.1                 
    ##  [57] reshape2_1.4.4              abind_1.4-5                
    ##  [59] httpuv_1.6.3                tools_4.1.0                
    ##  [61] ellipsis_0.3.2              spatstat.core_2.3-0        
    ##  [63] RColorBrewer_1.1-2          ggridges_0.5.3             
    ##  [65] plyr_1.8.6                  zlibbioc_1.40.0            
    ##  [67] purrr_0.3.4                 RCurl_1.98-1.5             
    ##  [69] rpart_4.1-15                deldir_1.0-6               
    ##  [71] pbapply_1.5-0               cowplot_1.1.1              
    ##  [73] zoo_1.8-9                   SummarizedExperiment_1.24.0
    ##  [75] ggrepel_0.9.1               cluster_2.1.2              
    ##  [77] magrittr_2.0.1              RSpectra_0.16-0            
    ##  [79] data.table_1.14.2           scattermore_0.7            
    ##  [81] lmtest_0.9-38               RANN_2.6.1                 
    ##  [83] SnowballC_0.7.0             fitdistrplus_1.1-6         
    ##  [85] matrixStats_0.61.0          hms_1.1.1                  
    ##  [87] mime_0.12                   evaluate_0.14              
    ##  [89] xtable_1.8-4                XML_3.99-0.8               
    ##  [91] sparsesvd_0.2               gridExtra_2.3              
    ##  [93] compiler_4.1.0              tibble_3.1.6               
    ##  [95] KernSmooth_2.23-20          crayon_1.4.2               
    ##  [97] R.oo_1.24.0                 htmltools_0.5.2            
    ##  [99] tzdb_0.2.0                  mgcv_1.8-38                
    ## [101] later_1.3.0                 tidyr_1.1.4                
    ## [103] DBI_1.1.1                   tweenr_1.0.2               
    ## [105] MASS_7.3-54                 readr_2.0.2                
    ## [107] Matrix_1.3-4                cli_3.3.0                  
    ## [109] R.methodsS3_1.8.1           parallel_4.1.0             
    ## [111] igraph_1.2.8                pkgconfig_2.0.3            
    ## [113] TFMPvalue_0.0.8             GenomicAlignments_1.30.0   
    ## [115] plotly_4.10.0               spatstat.sparse_2.0-0      
    ## [117] xml2_1.3.2                  annotate_1.72.0            
    ## [119] svglite_2.0.0               DirichletMultinomial_1.36.0
    ## [121] webshot_0.5.2               rvest_1.0.2                
    ## [123] stringr_1.4.0               digest_0.6.28              
    ## [125] pracma_2.3.3                sctransform_0.3.3          
    ## [127] RcppAnnoy_0.0.19            CNEr_1.30.0                
    ## [129] spatstat.data_2.1-0         rmarkdown_2.11             
    ## [131] leiden_0.3.9                fastmatch_1.1-3            
    ## [133] uwot_0.1.10                 restfulr_0.0.13            
    ## [135] gtools_3.9.2                shiny_1.7.1                
    ## [137] Rsamtools_2.10.0            rjson_0.2.20               
    ## [139] lifecycle_1.0.1             nlme_3.1-153               
    ## [141] jsonlite_1.7.2              viridisLite_0.4.0          
    ## [143] fansi_0.5.0                 pillar_1.6.4               
    ## [145] lattice_0.20-45             GO.db_3.14.0               
    ## [147] KEGGREST_1.34.0             fastmap_1.1.0              
    ## [149] httr_1.4.2                  survival_3.2-13            
    ## [151] glue_1.6.2                  qlcMatrix_0.9.7            
    ## [153] png_0.1-7                   bit_4.0.4                  
    ## [155] ggforce_0.3.3               stringi_1.7.5              
    ## [157] blob_1.2.2                  caTools_1.18.2             
    ## [159] memoise_2.0.0               irlba_2.3.3                
    ## [161] future.apply_1.8.1
