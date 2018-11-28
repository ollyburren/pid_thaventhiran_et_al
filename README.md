# pid_thaventhiran_et_al
Code for Primary Immune Deficiency (PID) analysis presented in  Thaventhiran et al.
This code is for illustrative purposes only and will not run without necessary resource files.

## R session info
> sessionInfo()
R version 3.3.3 (2017-03-06)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: Scientific Linux 7.5 (Nitrogen)

locale:
 [1] LC_CTYPE=en_GB.UTF-8          LC_NUMERIC=C
 [3] LC_TIME=en_GB.UTF-8           LC_COLLATE=en_GB.UTF-8
 [5] LC_MONETARY=en_GB.UTF-8       LC_MESSAGES=en_GB.UTF-8
 [7] LC_PAPER=en_GB.UTF-8          LC_NAME=en_GB.UTF-8
 [9] LC_ADDRESS=en_GB.UTF-8        LC_TELEPHONE=en_GB.UTF-8
[11] LC_MEASUREMENT=en_GB.UTF-8    LC_IDENTIFICATION=en_GB.UTF-8

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils
 [8] datasets  methods   base

other attached packages:
 [1] xlsx_0.5.7                 xlsxjars_0.6.1
 [3] rJava_0.9-9                wgsea_1.8
 [5] snpStats_1.24.0            Matrix_1.2-8
 [7] survival_2.41-3            simGWAS_0.1-0
 [9] rtracklayer_1.32.2         reshape2_1.4.3
[11] RColorBrewer_1.1-2         rCOGS_0.0.0.9000
[13] optparse_1.4.4             Gviz_1.18.2
[15] ggrepel_0.7.0              GenomicInteractions_1.8.1
[17] InteractionSet_1.2.1       SummarizedExperiment_1.4.0
[19] Biobase_2.34.0             GenomicRanges_1.26.4
[21] GenomeInfoDb_1.10.3        IRanges_2.8.2
[23] S4Vectors_0.12.2           BiocGenerics_0.20.0
[25] cowplot_0.9.2              ggplot2_2.2.1
[27] biomaRt_2.30.0             magrittr_1.5
[29] devtools_1.13.5            data.table_1.11.2

loaded via a namespace (and not attached):
 [1] bitops_1.0-6                  matrixStats_0.53.1
 [3] bit64_0.9-7                   httr_1.3.1
 [5] tools_3.3.3                   backports_1.1.2
 [7] R6_2.2.2                      rpart_4.1-10
 [9] Hmisc_4.1-1                   DBI_1.0.0
[11] lazyeval_0.2.1                colorspace_1.3-2
[13] nnet_7.3-12                   withr_2.1.1
[15] gridExtra_2.3                 bit_1.1-12
[17] htmlTable_1.11.2              scales_0.5.0
[19] checkmate_1.8.5               mvtnorm_1.0-7
[21] stringr_1.3.1                 digest_0.6.15
[23] Rsamtools_1.26.2              foreign_0.8-69
[25] XVector_0.14.1                base64enc_0.1-3
[27] dichromat_2.0-0               pkgconfig_2.0.1
[29] htmltools_0.3.6               ensembldb_1.6.2
[31] BSgenome_1.42.0               htmlwidgets_1.0
[33] rlang_0.2.0                   rstudioapi_0.7
[35] RSQLite_2.1.1                 BiocInstaller_1.24.0
[37] shiny_1.0.5                   bindr_0.1
[39] combinat_0.0-8                BiocParallel_1.8.2
[41] acepack_1.4.1                 dplyr_0.7.4
[43] VariantAnnotation_1.20.3      RCurl_1.95-4.10
[45] Formula_1.2-2                 Rcpp_0.12.16
[47] munsell_0.4.3                 stringi_1.2.4
[49] yaml_2.2.0                    zlibbioc_1.20.0
[51] plyr_1.8.4                    AnnotationHub_2.6.5
[53] blob_1.1.1                    lattice_0.20-34
[55] Biostrings_2.42.1             splines_3.3.3
[57] GenomicFeatures_1.26.4        knitr_1.20
[59] pillar_1.2.1                  igraph_1.1.2
[61] corpcor_1.6.9                 XML_3.98-1.11
[63] glue_1.2.0                    biovizBase_1.22.0
[65] latticeExtra_0.6-28           httpuv_1.3.6.2
[67] getopt_1.20.2                 gtable_0.2.0
[69] assertthat_0.2.0              mime_0.5
[71] xtable_1.8-2                  tibble_1.4.2
[73] GenomicAlignments_1.10.1      AnnotationDbi_1.36.2
[75] memoise_1.1.0                 bindrcpp_0.2
[77] cluster_2.0.6                 interactiveDisplayBase_1.12.0
