# ============================================================
# 03_make_figures.R
# Creates figures from saved DESeq2 objects:
# - PCA plot (VST)
# - Volcano plot
# Saves to results/figures/
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(readr)
  library(dplyr)
})

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

make_figures <- function(prefix) {
  message("\n=== Figures: ", prefix, " ===")
  
  dds_path <- file.path("results/tables", paste0(prefix, "_dds.rds"))
  res_path <- file.path("results/tables", paste0(prefix, "_deseq2_results.csv"))
  
  if (!file.exists(dds_path)) stop("Missing: ", dds_path)
  if (!file.exists(res_path)) stop("Missing: ", res_path)
  
  dds <- readRDS(dds_path)
  res <- read_csv(res_path, show_col_types = FALSE)
  
  # ----------------------------
  # PCA (VST)
  # ----------------------------
  vsd <- vst(dds, blind = TRUE)
  pca <- prcomp(t(assay(vsd)))
  
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    condition = colData(dds)$condition
  )
  
  p_pca <- ggplot(pca_df, aes(PC1, PC2, shape = condition)) +
    geom_point(size = 2.2) +
    theme_minimal() +
    labs(title = paste0(prefix, ": PCA (VST-transformed counts)"))
  
  pca_out <- file.path("results/figures", paste0(prefix, "_pca.png"))
  ggsave(pca_out, p_pca, width = 7, height = 5)
  message("Saved: ", pca_out)
  
  # ----------------------------
  # Volcano
  # ----------------------------
  res2 <- res %>%
    mutate(
      neglog10padj = -log10(padj),
      neglog10padj = ifelse(is.finite(neglog10padj), neglog10padj, NA_real_)
    )
  
  p_volcano <- ggplot(res2, aes(x = log2FoldChange, y = neglog10padj)) +
    geom_point(alpha = 0.35) +
    theme_minimal() +
    labs(
      title = paste0(prefix, ": Volcano plot"),
      x = "log2 fold change (tumor vs normal)",
      y = "-log10 adjusted p-value"
    )
  
  volcano_out <- file.path("results/figures", paste0(prefix, "_volcano.png"))
  ggsave(volcano_out, p_volcano, width = 7, height = 5)
  message("Saved: ", volcano_out)
  
  invisible(TRUE)
}

make_figures("tcga_brca")
make_figures("tcga_luad")

message("\nAll figures created.")
