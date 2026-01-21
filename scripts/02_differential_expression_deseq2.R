# ============================================================
# 02_differential_expression_deseq2.R
# Differential expression with DESeq2 (OOM-safe)
# - Keeps all normal samples
# - Subsamples tumor samples
# - Prefilters low-expression genes
# - Saves checkpoints
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(readr)
})

# ---- Thread limits (important for stability)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# DESeq2 runner
# ------------------------------------------------------------
run_deseq2 <- function(
    prefix,
    max_tumor = 300,
    min_count = 10,
    min_samples = 10
) {
  message("\n=== DESeq2: ", prefix, " ===")
  
  # ---- Load data
  counts <- readRDS(file.path("data/processed", paste0(prefix, "_counts.rds")))
  coldata <- readRDS(file.path("data/processed", paste0(prefix, "_coldata.rds")))
  
  # ensure alignment
  coldata <- coldata[colnames(counts), , drop = FALSE]
  
  message(
    "Samples: ", nrow(coldata),
    " | tumor=", sum(coldata$condition == "tumor"),
    " normal=", sum(coldata$condition == "normal")
  )
  
  # ---- Subsample tumors (RAM protection)
  set.seed(42)
  
  tumor_ids  <- rownames(coldata)[coldata$condition == "tumor"]
  normal_ids <- rownames(coldata)[coldata$condition == "normal"]
  
  if (length(tumor_ids) > max_tumor) {
    tumor_ids <- sample(tumor_ids, max_tumor)
  }
  
  keep_ids <- c(normal_ids, tumor_ids)
  
  counts_sub <- counts[, keep_ids, drop = FALSE]
  coldata_sub <- coldata[keep_ids, , drop = FALSE]
  
  message(
    "Samples used: ", nrow(coldata_sub),
    " | tumor=", sum(coldata_sub$condition == "tumor"),
    " normal=", sum(coldata_sub$condition == "normal")
  )
  
  # ---- Gene prefiltering (huge RAM/time saver)
  keep_genes <- rowSums(counts_sub >= min_count) >= min_samples
  counts_sub <- counts_sub[keep_genes, , drop = FALSE]
  
  message("Genes kept after prefilter: ", nrow(counts_sub))
  
  # ---- DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData   = coldata_sub,
    design    = ~ condition
  )
  
  dds <- DESeq(dds)
  
  # ---- Checkpoint save (important)
  saveRDS(dds, file.path("results/tables", paste0(prefix, "_dds.rds")))
  
  # ---- Results
  res <- results(dds, contrast = c("condition", "tumor", "normal"))
  res_df <- as.data.frame(res) |>
    rownames_to_column("gene_id") |>
    arrange(padj)
  
  write_csv(
    res_df,
    file.path("results/tables", paste0(prefix, "_deseq2_results.csv"))
  )
  
  message(
    "Significant genes (padj < 0.05): ",
    sum(res_df$padj < 0.05, na.rm = TRUE)
  )
  
  invisible(res_df)
}

# ------------------------------------------------------------
# Run analyses
# ------------------------------------------------------------
brca_res <- run_deseq2("tcga_brca", max_tumor = 300)
luad_res <- run_deseq2("tcga_luad", max_tumor = 300)

message("\n=== DESeq2 analysis completed successfully ===")
