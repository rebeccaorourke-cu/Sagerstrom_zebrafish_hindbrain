Cluster Naming and plots neural HB13hpf R Notebook
================

# 1. libraries and palette

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

    ## Warning: package 'ggplot2' was built under R version 4.1.2

    ## Warning: package 'limma' was built under R version 4.1.3

    ## Warning: package 'patchwork' was built under R version 4.1.2

    ## Warning: package 'TFBSTools' was built under R version 4.1.1

    ## Warning: package 'motifmatchr' was built under R version 4.1.1

``` r
options(future.globals.maxSize = 4000 * 1024^2)
```

``` r
mypal <- pal_igv(palette = "default",alpha = 1)(35)
```

# 2. Read data

``` r
seurat <- readRDS(file = "RDSfiles/HB13hpf_neural.RDS")
DefaultAssay(seurat) <- "SCT"
Idents(seurat) <- "wsnn_res.8"
DimPlot(seurat, reduction = "wnn.umap", label = T, repel = T) + scale_color_igv()
```

![](Figure1and5_files/figure-gfm/readData-1.png)<!-- -->

# 3. Rename Idents

``` r
Idents(seurat) <- "wsnn_res.8"
seurat <- RenameIdents(seurat,
                       "0" = "CHB.1",
                       "1" = "r3",
                       "2" = "CHB.2",
                       "3" = "MHB.1",
                       "4" = "SC.1",
                       "5" = "r5",
                       "6" = "r1",
                       "7" = "CHB.3",
                       "8" = "r4",
                       "9" = "r6",
                       "10" = "MB.1",
                       "11" = "Neuron",
                       "12" = "r2",
                       "13" = "FB.1",
                       "14" = "MB.2",
                       "15" = "MHB.2",
                       "16" = "MHB.3",
                       "17" = "MHB.4",
                       "18" = "SC.2",
                       "19" = "MB.3",
                       "20" = "FB.2",
                       "21" = "r1 & r2",
                       "22" = "Ciliated",
                       "23" = "FB.3",
                       "24" = "FB.4",
                       "25" = "MHB.5",
                       "26" = "SC.3")
levels(seurat) <- c("FB.1","FB.2","FB.3","FB.4","MB.1","MB.2","MB.3","MHB.1","MHB.2","MHB.3",
                    "MHB.4","MHB.5","r1","r2","r1 & r2","r3","r4","r5","r6","CHB.1","CHB.2",
                    "CHB.3","SC.1","SC.2","SC.3","Neuron","Ciliated")
umapPlot <- DimPlot(seurat, reduction = "wnn.umap", pt.size = 3) + scale_color_igv() + guides(color = guide_legend(override.aes = list(size=4), ncol=1) )
umapPlot
```

![](Figure1and5_files/figure-gfm/RenameIdents-1.png)<!-- -->

``` r
seurat$Clusters <- Idents(seurat)
#ggsave(filename = "../results/Fig1_HB13hpf_umapPlot.png", plot = umapPlot)
```

``` r
saveRDS(seurat, file = "RDSfiles/HB13hpf_neural.RDS")
```

# 4. Find DE genes

``` r
All.markers <- FindAllMarkers(seurat, only.pos = T, verbose = F)
write.table(All.markers, file = "../results/DataS2_Fig1_HB13hpf_markers.txt", sep = "\t", quote = F, col.names = NA)
```

``` r
top5.pval <- All.markers %>% group_by(cluster) %>% top_n(n=-5, wt = p_val) %>% top_n(n=5, wt = avg_log2FC)
top5.pval
```

    ## # A tibble: 135 × 7
    ## # Groups:   cluster [27]
    ##        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene           
    ##        <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>          
    ##  1 2.48e- 66      0.800 0.459 0.008 4.54e- 62 FB.1    rx1            
    ##  2 3.31e- 58      1.38  0.649 0.033 6.06e- 54 FB.1    aldh1a3        
    ##  3 4.05e- 45      1.10  0.459 0.02  7.42e- 41 FB.1    rorb           
    ##  4 5.97e- 44      0.976 0.568 0.035 1.09e- 39 FB.1    six3b          
    ##  5 5.36e- 40      0.344 0.216 0.001 9.81e- 36 FB.1    hmx4           
    ##  6 6.25e-112      2.42  1     0.028 1.14e-107 FB.2    si:ch73-215f7.1
    ##  7 6.57e- 89      1.21  0.69  0.013 1.20e- 84 FB.2    foxg1a         
    ##  8 6.22e- 65      1.66  0.724 0.029 1.14e- 60 FB.2    emx3           
    ##  9 2.99e- 47      0.554 0.379 0.008 5.48e- 43 FB.2    dmrta2         
    ## 10 2.01e- 42      0.730 0.379 0.01  3.68e- 38 FB.2    dlx3b          
    ## # … with 125 more rows

# 5. Plots

## 5.1 dotplot

``` r
dotPlot <- DotPlot(seurat, features = unique(top5.pval$gene)) + RotatedAxis()
dotPlot
```

![](Figure1and5_files/figure-gfm/dotplot-1.png)<!-- -->

``` r
#ggsave(filename = "../results/Fig1_HB13hpf_dotPlot.png", plot = dotPlot)
```

## 5.2 heatmap

``` r
heatmapPlot <- DoHeatmap(seurat, features = unique(top5.pval$gene), group.colors = mypal, 
                         size = 5, angle = 45) + 
  guides(color = FALSE) +
  theme(axis.text = element_blank())
```

    ## Warning: The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as
    ## of ggplot2 3.3.4.

``` r
heatmapPlot
```

![](Figure1and5_files/figure-gfm/heatmap-1.png)<!-- -->

``` r
ggsave(filename = "../results/Fig5_HB13hpf_heatmapPlot.png", plot = heatmapPlot) ## heatmap included in Fig 5 not Fig 1
```

    ## Saving 7 x 7 in image

# Build Cluster Tree

``` r
seurat <- BuildClusterTree(seurat, graph = "wsnn")
```

``` r
p <- PlotClusterTree(seurat,direction = "rightwards")
```

![](Figure1and5_files/figure-gfm/plotclustertree-1.png)<!-- -->

``` r
p
```

    ## NULL

``` r
data.tree <- Tool(object = seurat, slot = "BuildClusterTree")
pdf(file = "../results/Fig5_HB13hpf_clustertree_unrooted_vrs1.pdf")
p1 <- ape::plot.phylo(x = data.tree,  use.edge.length = FALSE, cex = 0.8, type = "u", 
                      align.tip.label = TRUE, label.offset = 0.5, lab4ut = "axial")
p1
```

    ## $type
    ## [1] "unrooted"
    ## 
    ## $use.edge.length
    ## [1] FALSE
    ## 
    ## $node.pos
    ## NULL
    ## 
    ## $node.depth
    ## [1] 1
    ## 
    ## $show.tip.label
    ## [1] TRUE
    ## 
    ## $show.node.label
    ## [1] FALSE
    ## 
    ## $font
    ##  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
    ## 
    ## $cex
    ##  [1] 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8
    ## [20] 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8
    ## 
    ## $adj
    ##  [1] 0 0 0 0 0 1 1 1 1 1 1 1 0 1 1 1 1 1 0 0 0 0 1 0 1 0 0
    ## 
    ## $srt
    ## [1] 0
    ## 
    ## $no.margin
    ## [1] FALSE
    ## 
    ## $label.offset
    ## [1] 0.5
    ## 
    ## $x.lim
    ## [1] -1.113558  8.619067
    ## 
    ## $y.lim
    ## [1] -1.113558 10.779861
    ## 
    ## $direction
    ## [1] "rightwards"
    ## 
    ## $tip.color
    ##  [1] "black" "black" "black" "black" "black" "black" "black" "black" "black"
    ## [10] "black" "black" "black" "black" "black" "black" "black" "black" "black"
    ## [19] "black" "black" "black" "black" "black" "black" "black" "black" "black"
    ## 
    ## $Ntip
    ## [1] 27
    ## 
    ## $Nnode
    ## [1] 26
    ## 
    ## $root.time
    ## NULL
    ## 
    ## $align.tip.label
    ## [1] FALSE

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
pdf(file = "../results/Fig5_HB13hpf_clustertree_unrooted_vrs2.pdf")
p1 <- ape::plot.phylo(x = data.tree,  use.edge.length = FALSE, cex = 0.8, type = "u", 
                      align.tip.label = TRUE, label.offset = 0.5)
p1
```

    ## $type
    ## [1] "unrooted"
    ## 
    ## $use.edge.length
    ## [1] FALSE
    ## 
    ## $node.pos
    ## NULL
    ## 
    ## $node.depth
    ## [1] 1
    ## 
    ## $show.tip.label
    ## [1] TRUE
    ## 
    ## $show.node.label
    ## [1] FALSE
    ## 
    ## $font
    ## [1] 3
    ## 
    ## $cex
    ## [1] 0.8
    ## 
    ## $adj
    ## [1] 0
    ## 
    ## $srt
    ## [1] 0
    ## 
    ## $no.margin
    ## [1] FALSE
    ## 
    ## $label.offset
    ## [1] 0.5
    ## 
    ## $x.lim
    ## [1] -1.113558  8.619067
    ## 
    ## $y.lim
    ## [1] -1.113558 10.779861
    ## 
    ## $direction
    ## [1] "rightwards"
    ## 
    ## $tip.color
    ## [1] "black"
    ## 
    ## $Ntip
    ## [1] 27
    ## 
    ## $Nnode
    ## [1] 26
    ## 
    ## $root.time
    ## NULL
    ## 
    ## $align.tip.label
    ## [1] FALSE

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## 5.3 combine umap, dotplot and heatmap

``` r
# combined <- (((umapPlot + 
#       theme(legend.text = element_text(size = 20))) |
#      (heatmapPlot +
#         theme(axis.text = element_blank()))) / 
#     (dotPlot + 
#        theme(axis.text.x = element_text(size = 10),
#              axis.text.y = element_text(size = 12),
#              legend.position = "bottom"))) +
#     plot_layout(heights = c(2.1,1))
# combined
# ggsave(filename = "../results/Fig1_HB13hpf_combined2Plot.png", plot = combined)
```

``` r
#genes1 <- FeaturePlot(seurat, features = c("sox3","neurod4"), reduction = "wnn.umap", keep.scale = "all")  
#genes1
#ggsave(filename = "../results/Fig1_HB13hpf_genes1_scaledToMax.png", plot = genes1)
```

``` r
genes1 <- FeaturePlot(seurat, features = c("sox3","neurod4"), reduction = "wnn.umap", max.cutoff = 1)  
genes1
```

![](Figure1and5_files/figure-gfm/featureplot-1.png)<!-- -->

``` r
#ggsave(filename = "../results/Fig1_HB13hpf_genes1_scaledToMin.png", plot = genes1)
```

``` r
genelist <- c("egr2b","mafba","sema3ab","fgf3","cyp26b1","vgll3","tgfbr2b","irx1b")
```

``` r
genes2 <- FeaturePlot(seurat, features = genelist, reduction = "wnn.umap", max.cutoff = 1, combine = FALSE)
violinPlot <- VlnPlot(seurat, features = genelist, same.y.lims = TRUE, combine = FALSE)
```

``` r
genePlusViolin <- genes2[[1]] + violinPlot[[1]] + NoLegend() + genes2[[2]] + violinPlot[[2]] + NoLegend() + 
  genes2[[3]] + violinPlot[[3]] + NoLegend() + genes2[[4]] + violinPlot[[4]] + NoLegend() +
  genes2[[5]] + violinPlot[[5]] + NoLegend() + genes2[[6]] + violinPlot[[6]] + NoLegend() +
  genes2[[7]] + violinPlot[[7]] + NoLegend() + genes2[[8]] + violinPlot[[8]] + NoLegend() +
  plot_layout(ncol = 4, widths = c(1,2,1,2)) #& scale_color_gradient(low = "gray", high = "blue", limits = c(0,1))
genePlusViolin
```

![](Figure1and5_files/figure-gfm/geneplusvln-1.png)<!-- -->

``` r
#ggsave(filename = "../results/Fig1_HB13hpf_genePlusViolin.png", plot = genePlusViolin)
```

``` r
genes3 <- FeaturePlot(seurat, features = c("casz1","zic2a","ntn1a","sp8b"), reduction = "wnn.umap", ncol = 4, max.cutoff = 1) #&
  #scale_color_gradient(low = "gray", high = "blue", limits = c(0,2.4)) 
genes3
```

![](Figure1and5_files/figure-gfm/featureplot3-1.png)<!-- -->

``` r
combined <- 
  (((umapPlot +
      theme(legend.text = element_text(size = 20))) /
      plot_spacer() /
     (genes1 + plot_layout(guides = "collect")) /
     plot_spacer() + 
     plot_layout(heights = c(2,0.1,0.7,3.1))
  ) |
  ((dotPlot +
      theme(axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 12),
            legend.position = "bottom")) /
     plot_spacer() /
     (genePlusViolin + plot_layout(guides = "collect")) /
     plot_spacer() /
     (genes3 + plot_layout(guides = "collect")) +
     plot_layout(heights = c(1,0.1,4,0.1,0.7))
  )) + 
  plot_layout(widths = c(1,2.2))
combined
```

![](Figure1and5_files/figure-gfm/combined2-1.png)<!-- -->

``` r
ggsave(filename = "../results/Fig1_HB13hpf_combinedPlot.png", plot = combined)
```

    ## Saving 25 x 25 in image

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
    ##  [5] patchwork_1.1.2                     JASPAR2020_0.99.10                 
    ##  [7] limma_3.50.3                        ggsci_2.9                          
    ##  [9] ggplot2_3.4.0                       dplyr_1.0.7                        
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
    ##  [13] ragg_1.2.4                  codetools_0.2-18           
    ##  [15] ica_1.0-2                   future_1.26.1              
    ##  [17] miniUI_0.1.1.1              withr_2.5.0                
    ##  [19] colorspace_2.0-2            highr_0.9                  
    ##  [21] rstudioapi_0.13             ROCR_1.0-11                
    ##  [23] tensor_1.5                  listenv_0.8.0              
    ##  [25] labeling_0.4.2              MatrixGenerics_1.6.0       
    ##  [27] slam_0.1-48                 GenomeInfoDbData_1.2.7     
    ##  [29] polyclip_1.10-0             bit64_4.0.5                
    ##  [31] farver_2.1.0                parallelly_1.32.0          
    ##  [33] vctrs_0.5.0                 generics_0.1.1             
    ##  [35] xfun_0.27                   lsa_0.73.2                 
    ##  [37] ggseqlogo_0.1               R6_2.5.1                   
    ##  [39] bitops_1.0-7                spatstat.utils_2.2-0       
    ##  [41] cachem_1.0.6                DelayedArray_0.20.0        
    ##  [43] assertthat_0.2.1            promises_1.2.0.1           
    ##  [45] BiocIO_1.4.0                scales_1.2.1               
    ##  [47] gtable_0.3.0                globals_0.15.1             
    ##  [49] goftest_1.2-3               seqLogo_1.60.0             
    ##  [51] rlang_1.0.6                 systemfonts_1.0.4          
    ##  [53] RcppRoll_0.3.0              splines_4.1.0              
    ##  [55] lazyeval_0.2.2              spatstat.geom_2.3-0        
    ##  [57] yaml_2.2.1                  reshape2_1.4.4             
    ##  [59] abind_1.4-5                 httpuv_1.6.3               
    ##  [61] tools_4.1.0                 ellipsis_0.3.2             
    ##  [63] spatstat.core_2.3-0         RColorBrewer_1.1-2         
    ##  [65] ggridges_0.5.3              plyr_1.8.6                 
    ##  [67] zlibbioc_1.40.0             purrr_0.3.4                
    ##  [69] RCurl_1.98-1.5              rpart_4.1-15               
    ##  [71] deldir_1.0-6                pbapply_1.5-0              
    ##  [73] cowplot_1.1.1               zoo_1.8-9                  
    ##  [75] SummarizedExperiment_1.24.0 ggrepel_0.9.1              
    ##  [77] cluster_2.1.2               magrittr_2.0.1             
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
    ## [107] Matrix_1.3-4                cli_3.4.1                  
    ## [109] R.methodsS3_1.8.1           parallel_4.1.0             
    ## [111] igraph_1.2.8                pkgconfig_2.0.3            
    ## [113] TFMPvalue_0.0.8             GenomicAlignments_1.30.0   
    ## [115] plotly_4.10.0               spatstat.sparse_2.0-0      
    ## [117] xml2_1.3.3                  annotate_1.72.0            
    ## [119] svglite_2.1.0               DirichletMultinomial_1.36.0
    ## [121] webshot_0.5.4               rvest_1.0.3                
    ## [123] stringr_1.4.0               digest_0.6.28              
    ## [125] pracma_2.3.3                sctransform_0.3.3          
    ## [127] RcppAnnoy_0.0.19            CNEr_1.30.0                
    ## [129] spatstat.data_2.1-0         rmarkdown_2.11             
    ## [131] leiden_0.3.9                fastmatch_1.1-3            
    ## [133] uwot_0.1.10                 restfulr_0.0.13            
    ## [135] gtools_3.9.2                shiny_1.7.1                
    ## [137] Rsamtools_2.10.0            rjson_0.2.20               
    ## [139] lifecycle_1.0.3             nlme_3.1-153               
    ## [141] jsonlite_1.7.2              viridisLite_0.4.0          
    ## [143] fansi_0.5.0                 pillar_1.6.4               
    ## [145] lattice_0.20-45             GO.db_3.14.0               
    ## [147] KEGGREST_1.34.0             fastmap_1.1.0              
    ## [149] httr_1.4.2                  survival_3.2-13            
    ## [151] glue_1.6.2                  qlcMatrix_0.9.7            
    ## [153] png_0.1-7                   bit_4.0.4                  
    ## [155] ggforce_0.3.3               stringi_1.7.5              
    ## [157] blob_1.2.2                  textshaping_0.3.6          
    ## [159] caTools_1.18.2              memoise_2.0.0              
    ## [161] ape_5.6-2                   irlba_2.3.3                
    ## [163] future.apply_1.8.1
