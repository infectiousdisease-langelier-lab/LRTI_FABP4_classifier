This is the repository for testing the gene *FABP4* as an age-agnostic inverse-biomarker for distinguishing lower respiratory tract infection (LRTI) from non-LRTI.

## Methods

### Data preprocessing

The gene count and metadata files of the pediatric samples were downloaded from NCBI GEO accession number [GSE212532](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212532) (1). The gene count and metadata files are available in the folder [adult_count](adult_count) (2). Alignment was performed using kallisto (3) against the protein-coding and lncRNA transcripts from the human referece genome (Ensembl 99 release). The adult samples have been filtered to remove samples with proportion of rRNA lower than 80%.

### Data analysis

All analyses were performed in R (v4.3.2), and can be reproduced with the code [FABP4_classifier.R](FABP4_classifier.R).

For DE analysis, we compared "Definite" (those with LRTI) against "No Evidence" (those without LRTI) patients. First, we filtered for genes with at least 10 counts in at least 20% of the samples. We then performed DE testing with the limma-voom procedure (with `normalize.method="quantile"`) from the limma package (v3.58.1) (4,5). We did not take into account any covariates in the pediatric analysis, but we did control for batch number in the adult analysis. P-values were adjusted with Benjamini-Hochberg correction. The DE results for the pediatric and adult cohorts are available in [DE_children.csv](DE_children.csv) and [DE_adults.csv](DE_adults.csv), respectively.

To test the performance of *FABP4* as a biomarker for LRTI, we applied five-fold cross-validation separately on the pediatric samples and the adult samples. For each cohort, we randomly split the samples into five folds, such that the number of "Definite" and "No Evidence" samples are roughly equal across all 5 folds. Next, for each test fold, we filtered for genes with at least 10 counts in at least 20% of the training folds' samples, and normalized the gene counts with `varianceStabilizingTransformation()` (DESeq2 package v1.42.0) (6). (More details on the normalization step is provided [below](#normalization-of-fabp4-gene-expression).) For each test fold, the performance of *FABP4* was calculated using the function `roc()` from the pROC package (v1.18.5) (7). The reported area under the receiver operating characteristic curve was calculated from the mean and standard deviation of the 5 corresponding area under the curve values. The mean ROC curve was calculated from five linearly interpolated ROC curves, one for each test fold.

### Normalization of *FABP4* gene expression

For each of the 5 folds, *FABP4* gene expression level was normalized using the DESeq2 package. First, we filtered for genes with at least 10 counts in at least 20% of the training samples. Next, we calculated the training samples' size factors and dispersions with the function `estimateSizeFactors()` and `estimateDispersions()`, respectively. Third, we normalized the training samples' gene expression with the function `varianceStabilizingTransformation()`, rounded to two decimal places. For the test samples, we filtered for genes that were analyzed in the training samples. Next, we calculated their size factors, but set their dispersions to be equal to the training samples' dispersions. Then, we normalized the test samples' gene expression with `varianceStabilizingTransformation(blind=FALSE)` (so that the function used the preset dispersions), and rounded to two decimal places. The normalized *FABP4* gene expression is directly used as input to the pROC package.

## Code

FABP4_classifier.R: code used to analyze data and produce the figures.
* Figure 1A: `volcano_children.png`
* Figure 1B: `mean_roc.svg`

## SessionInfo

```
R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.6.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: US/Pacific
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] pROC_1.18.5                 DESeq2_1.42.0               SummarizedExperiment_1.32.0 Biobase_2.62.0              MatrixGenerics_1.14.0      
 [6] matrixStats_1.2.0           GenomicRanges_1.54.1        GenomeInfoDb_1.38.5         IRanges_2.36.0              S4Vectors_0.40.2           
[11] BiocGenerics_0.48.1         limma_3.58.1                patchwork_1.2.0             lubridate_1.9.3             forcats_1.0.0              
[16] stringr_1.5.1               dplyr_1.1.4                 purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                
[21] tibble_3.2.1                ggplot2_3.4.4               tidyverse_2.0.0            

loaded via a namespace (and not attached):
 [1] gtable_0.3.4            lattice_0.22-5          tzdb_0.4.0              vctrs_0.6.5             tools_4.3.2            
 [6] bitops_1.0-7            generics_0.1.3          parallel_4.3.2          fansi_1.0.6             pkgconfig_2.0.3        
[11] Matrix_1.6-5            lifecycle_1.0.4         GenomeInfoDbData_1.2.11 farver_2.1.1            compiler_4.3.2         
[16] textshaping_0.3.7       statmod_1.5.0           munsell_0.5.0           codetools_0.2-19        RCurl_1.98-1.14        
[21] pillar_1.9.0            crayon_1.5.2            BiocParallel_1.36.0     DelayedArray_0.28.0     abind_1.4-5            
[26] locfit_1.5-9.8          tidyselect_1.2.0        stringi_1.8.3           labeling_0.4.3          grid_4.3.2             
[31] colorspace_2.1-0        cli_3.6.2               SparseArray_1.2.3       magrittr_2.0.3          S4Arrays_1.2.0         
[36] utf8_1.2.4              withr_3.0.0             scales_1.3.0            timechange_0.3.0        XVector_0.42.0         
[41] ragg_1.2.7              hms_1.1.3               rlang_1.1.3             Rcpp_1.0.12             glue_1.7.0             
[46] rstudioapi_0.15.0       R6_2.5.1                plyr_1.8.9              systemfonts_1.0.5       zlibbioc_1.48.0        
```

## References

1. Mick, Eran, et al. "Integrated host/microbe metagenomics enables accurate lower respiratory tract infection diagnosis in critically ill children." Journal of Clinical Investigation 133.7 (2023): e165904.
2. Langelier, Charles, et al. "Integrating host response and unbiased microbe detection for lower respiratory tract infection diagnosis in critically ill adults." Proceedings of the National Academy of Sciences 115.52 (2018): E12353-E12362.
3. Bray, Nicolas L., et al. "Near-optimal probabilistic RNA-seq quantification." Nature biotechnology 34.5 (2016): 525-527.
4. Ritchie, Matthew E., et al. "limma powers differential expression analyses for RNA-sequencing and microarray studies." Nucleic acids research 43.7 (2015): e47-e47.
5. Law, Charity W., et al. "voom: Precision weights unlock linear model analysis tools for RNA-seq read counts." Genome biology 15.2 (2014): 1-17.
6. Love, Michael I., Wolfgang Huber, and Simon Anders. "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome biology 15.12 (2014): 1-21.
7. Robin, Xavier, et al. "pROC: an open-source package for R and S+ to analyze and compare ROC curves." BMC bioinformatics 12.1 (2011): 1-8.