# Reproducibility Guide

This document explains how to reproduce the full analysis from scratch.

The pipeline is designed to be lightweight, reproducible, and runnable on
a standard laptop.

---

## Requirements

### Software
- R (≥ 4.3)
- R packages:
  - recount3
  - SummarizedExperiment
  - DESeq2
  - tidyverse (dplyr, ggplot2, readr, tibble)

---

## Install Required R Packages

Start R or RStudio and run:

```r
install.packages("BiocManager")

BiocManager::install(c(
  "recount3",
  "SummarizedExperiment",
  "DESeq2"
))

install.packages(c(
  "dplyr",
  "ggplot2",
  "readr",
  "tibble"
))
```
However you can also install the packages using renv:
```r
renv::restore()
```

## Running the Analysis

### Step 1: Download and prepare data

```bash
Rscript scripts/01_download_prepare_recount3.R
```
This will:

- Query TCGA metadata via `recount3`
- Download required data on demand
- Create processed count matrices and sample metadata

Output:

```bash
data/processed/
├── tcga_brca_counts.rds
├── tcga_brca_coldata.rds
├── tcga_luad_counts.rds
├── tcga_luad_coldata.rds
```

### Step 2: Differential expression analysis

```bash
Rscript scripts/02_differential_expression_deseq2.R
```

This step:

- Retains all normal samples
- Subsamples tumor samples (for memory safety)
- Prefilters low-expression genes
- Runs DESeq2 and saves results

Output:

```bash
results/tables/
├── tcga_brca_deseq2_results.csv
├── tcga_luad_deseq2_results.csv
```

### Step 3: Generate figures

```bash
Rscript scripts/03_make_figures.R
```

This will generate:

- PCA plots (VST-transformed expression)
- Volcano plots (tumor vs normal)

Output:

```bash
results/figures/
├── tcga_brca_pca.png
├── tcga_brca_volcano.png
├── tcga_luad_pca.png
├── tcga_luad_volcano.png
```


