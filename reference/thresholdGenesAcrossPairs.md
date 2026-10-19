# Function that runs thresholding of genes by sample pairs

Function that runs thresholding of genes by sample pairs

## Usage

``` r
thresholdGenesAcrossPairs(
  seuratObject,
  uniqueSampleID = "sample",
  pairingColumn,
  percentThreshold = 0.01,
  completePairsOnly = FALSE,
  assay = "RNA"
)
```

## Arguments

- seuratObject:

  Bayesian Meta-Analysis object with basic bayesian meta-analysis run

- uniqueSampleID:

  Column in seurat metadata object of the unique sample ID

- pairingColumn:

  Column in seurat metadata object that samples are paired on

- percentThreshold:

  Percentage of cells that should express the gene across each pair to
  be included

- completePairsOnly:

  BOOL where TRUE subsets seuratObject to samples from only the complete
  pairs

- assay:

  assay to use

## Value

Vector of genes that pass thresholds
