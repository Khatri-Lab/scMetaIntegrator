# Function that runs detects complete pairs

Function that runs detects complete pairs

## Usage

``` r
confirmPairs(seuratObject, uniqueSampleID = "sample", pairingColumn)
```

## Arguments

- seuratObject:

  Seurat v5 Object

- uniqueSampleID:

  Column in seurat metadata object of the unique sample ID

- pairingColumn:

  Column in seurat metadata object that samples are paired on

## Value

boolean TRUE if two samples per pair present, FALSE if not
