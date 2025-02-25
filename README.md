# scMetaIntegrator
scMetaIntegrator: a meta-analysis approach to paired single-cell differential expression analysis\
Preprint: link will be available soon\

## Goal
scMetaIntegrator aims to provide a method for differential expresion analysis that can handle paired scRNA-seq data. 

## Summary
As single-cell RNA sequencing (scRNA-seq) studies become increasingly complex, traditional differential gene expression methods can fall short, particularly in paired repeated measures and matched cohort designs. Many existing approaches consider cells as independent samples, leading to false positives while ignoring inherent sampling structures. While pseudobulk methods can address this, they ignore intra-sample gene expression variability, resulting in high false negatives. We propose a novel meta-analysis approach that accounts for biological replicates and cell variability in paired scRNA-seq data. Our method, single-cell MetaIntegrator, yields robust effect size estimates and reproducible p-values, demonstrating effectiveness on both real and synthetic datasets. It also builds on our established MetaIntegrator package and its suite of visualization and analysis techniques.


## Installation
scMetaIntegrator is implemented in R. To install:
```
devtools::install_github("Khatri-Lab/scMetaIntegrator", ref = "main")
```

## Dependencies
scMetaIntegrator requires installation of Seurat version 5.0.0 or higher. It also requires the MetaIntegrator package which currently must be installed from source.
```
install.packages("https://cran.r-project.org/src/contrib/Archive/COCONUT/COCONUT_1.0.2.tar.gz",
                 repos=NULL, method="libcurl")
install.packages("https://cran.r-project.org/src/contrib/Archive/MetaIntegrator/MetaIntegrator_2.1.3.tar.gz",
                 repos=NULL, method="libcurl")
```
We are working on getting MetaIntegrator back on CRAN and building scMetaIntegrator into the existing package. Updates will follow soon.

## Vignettes
Vignettes are available in the vignettes directory of this repo.