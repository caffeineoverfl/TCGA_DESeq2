# Pan-Cancer Differential Gene Expression Analysis with TCGA RNA-seq Data

This project performs a reproducible differential gene expression analysis of
tumor versus normal samples using RNA-seq data from The Cancer Genome Atlas (TCGA).
Two cancer types are analyzed as case studies:

- **Breast Invasive Carcinoma (BRCA)**
- **Lung Adenocarcinoma (LUAD)**

The analysis uses the `recount3` resource for uniformly processed RNA-seq data
and the `DESeq2` framework for statistical modeling.

---

## Biological Motivation

Cancer development is associated with large-scale transcriptional changes.
By comparing tumor samples to matched normal tissues, we can identify genes
and pathways that are consistently dysregulated in cancer.

This project aims to:
- Validate tumor–normal separation using unsupervised methods (PCA)
- Identify differentially expressed genes using a robust statistical model
- Compare transcriptional patterns across different cancer types

---

## Key Results

### PCA (VST-normalized expression)
- Clear separation between tumor and normal samples in both BRCA and LUAD
- Tumor samples show increased heterogeneity, consistent with cancer biology

### Differential Expression
- Thousands of genes significantly differentially expressed (adjusted p < 0.05)
- Strong and symmetric up- and down-regulation patterns visible in volcano plots
- Results are consistent across cancer types, demonstrating robustness

---

## Repository Structure
```bash
├── data/
│ ├── raw/ # intentionally empty (data streamed via recount3)
│ └── processed/ # processed count matrices and metadata
├── scripts/
│ ├── 01_download_prepare_recount3.R
│ ├── 02_differential_expression_deseq2.R
│ └── 03_make_figures.R
├── results/
│ ├── tables/
│ │ └── *_deseq2_results.csv
│ └── figures/
│ ├── *_pca.png
│ └── *_volcano.png
├── README.md
└── Howtoreproduce.md
```

---

## Methods

- RNA-seq data accessed via `recount3`
- Count normalization and dispersion estimation using `DESeq2`
- Tumor samples subsampled for computational feasibility
- Lowly expressed genes filtered prior to modeling
- Visualization via PCA (VST) and volcano plots

---

## License

This project uses only publicly available data and is intended for educational
and research purposes.

