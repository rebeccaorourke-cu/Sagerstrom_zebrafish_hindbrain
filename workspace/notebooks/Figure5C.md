HB16hpf neural Cluster Tree R Notebook
================

# 1. libraries

``` r
library(Seurat)
```

    ## Attaching SeuratObject

``` r
library(Signac)
```

    ## 
    ## Attaching package: 'Signac'

    ## The following object is masked from 'package:Seurat':
    ## 
    ##     FoldChange

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
library(ggsci)
options(future.globals.maxSize = 4000 * 1024^2)
```

# 2. Read data

``` r
HB16hpf <- readRDS("RDSfiles/HB16hpf_neural.RDS")
Idents(HB16hpf) <- "Clusters"
DimPlot(HB16hpf, reduction = "wnn.umap") + scale_color_igv()
```

![](Figure5C_files/figure-gfm/readData-1.png)<!-- -->

# 3. Build Cluster Tree

``` r
HB16hpf <- BuildClusterTree(HB16hpf, graph = "wsnn")
```

# 4. Plot Cluster Tree

``` r
p <- PlotClusterTree(HB16hpf,direction = "rightwards")
```

![](Figure5C_files/figure-gfm/PlotTree-1.png)<!-- -->

``` r
p
```

    ## NULL

``` r
data.tree <- Tool(object = HB16hpf, slot = "BuildClusterTree")
png(filename = "../results/Fig5C_HB16hpf_clustertree_unrooted.png",width = 900, height = 900)
p1 <- ape::plot.phylo(x = data.tree,  use.edge.length = FALSE, type = "u")
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
    ## [1] 1
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
    ## [1] 0
    ## 
    ## $x.lim
    ## [1] -1.65425 10.58884
    ## 
    ## $y.lim
    ## [1] -1.65425 10.00905
    ## 
    ## $direction
    ## [1] "rightwards"
    ## 
    ## $tip.color
    ## [1] "black"
    ## 
    ## $Ntip
    ## [1] 23
    ## 
    ## $Nnode
    ## [1] 22
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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggsci_2.9          ggplot2_3.3.5      dplyr_1.0.7        Signac_1.2.1      
    ## [5] SeuratObject_4.0.4 Seurat_4.0.1      
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fastmatch_1.1-3        plyr_1.8.6             igraph_1.2.8          
    ##   [4] lazyeval_0.2.2         splines_4.1.0          BiocParallel_1.28.0   
    ##   [7] listenv_0.8.0          scattermore_0.7        SnowballC_0.7.0       
    ##  [10] GenomeInfoDb_1.30.0    digest_0.6.28          htmltools_0.5.2       
    ##  [13] fansi_0.5.0            magrittr_2.0.1         tensor_1.5            
    ##  [16] cluster_2.1.2          ROCR_1.0-11            globals_0.15.1        
    ##  [19] Biostrings_2.62.0      matrixStats_0.61.0     docopt_0.7.1          
    ##  [22] spatstat.sparse_2.0-0  colorspace_2.0-2       ggrepel_0.9.1         
    ##  [25] xfun_0.27              sparsesvd_0.2          crayon_1.4.2          
    ##  [28] RCurl_1.98-1.5         jsonlite_1.7.2         spatstat.data_2.1-0   
    ##  [31] ape_5.5                survival_3.2-13        zoo_1.8-9             
    ##  [34] glue_1.6.2             polyclip_1.10-0        gtable_0.3.0          
    ##  [37] zlibbioc_1.40.0        XVector_0.34.0         leiden_0.3.9          
    ##  [40] future.apply_1.8.1     BiocGenerics_0.40.0    abind_1.4-5           
    ##  [43] scales_1.1.1           DBI_1.1.1              miniUI_0.1.1.1        
    ##  [46] Rcpp_1.0.7             viridisLite_0.4.0      xtable_1.8-4          
    ##  [49] reticulate_1.22        spatstat.core_2.3-0    stats4_4.1.0          
    ##  [52] htmlwidgets_1.5.4      httr_1.4.2             RColorBrewer_1.1-2    
    ##  [55] ellipsis_0.3.2         ica_1.0-2              pkgconfig_2.0.3       
    ##  [58] farver_2.1.0           ggseqlogo_0.1          uwot_0.1.10           
    ##  [61] deldir_1.0-6           utf8_1.2.2             labeling_0.4.2        
    ##  [64] tidyselect_1.1.1       rlang_1.0.3            reshape2_1.4.4        
    ##  [67] later_1.3.0            munsell_0.5.0          tools_4.1.0           
    ##  [70] cli_3.3.0              generics_0.1.1         ggridges_0.5.3        
    ##  [73] evaluate_0.14          stringr_1.4.0          fastmap_1.1.0         
    ##  [76] yaml_2.2.1             goftest_1.2-3          knitr_1.36            
    ##  [79] fitdistrplus_1.1-6     purrr_0.3.4            RANN_2.6.1            
    ##  [82] pbapply_1.5-0          future_1.26.1          nlme_3.1-153          
    ##  [85] mime_0.12              slam_0.1-48            RcppRoll_0.3.0        
    ##  [88] compiler_4.1.0         rstudioapi_0.13        plotly_4.10.0         
    ##  [91] png_0.1-7              spatstat.utils_2.2-0   tibble_3.1.6          
    ##  [94] tweenr_1.0.2           stringi_1.7.5          highr_0.9             
    ##  [97] lattice_0.20-45        Matrix_1.3-4           vctrs_0.4.1           
    ## [100] pillar_1.6.4           lifecycle_1.0.1        spatstat.geom_2.3-0   
    ## [103] lmtest_0.9-38          RcppAnnoy_0.0.19       data.table_1.14.2     
    ## [106] cowplot_1.1.1          bitops_1.0-7           irlba_2.3.3           
    ## [109] httpuv_1.6.3           patchwork_1.1.1        GenomicRanges_1.46.0  
    ## [112] R6_2.5.1               promises_1.2.0.1       KernSmooth_2.23-20    
    ## [115] gridExtra_2.3          lsa_0.73.2             IRanges_2.28.0        
    ## [118] parallelly_1.32.0      codetools_0.2-18       MASS_7.3-54           
    ## [121] assertthat_0.2.1       withr_2.4.2            qlcMatrix_0.9.7       
    ## [124] sctransform_0.3.3      Rsamtools_2.10.0       S4Vectors_0.32.4      
    ## [127] GenomeInfoDbData_1.2.7 mgcv_1.8-38            parallel_4.1.0        
    ## [130] grid_4.1.0             rpart_4.1-15           tidyr_1.1.4           
    ## [133] rmarkdown_2.11         Rtsne_0.15             ggforce_0.3.3         
    ## [136] shiny_1.7.1
