# 01_download_prepare_recount3.R
suppressPackageStartupMessages({
  library(recount3)
  library(SummarizedExperiment)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
})

dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

find_col <- function(nms, patterns) {
  hits <- nms[str_detect(tolower(nms), paste0(patterns, collapse = "|"))]
  if (length(hits) == 0) NA_character_ else hits[1]
}

prepare_tcga_project <- function(project_id, out_prefix) {
  message("\n=== ", project_id, " ===")
  
  ap <- available_projects()
  proj <- ap[ap$project == project_id, ]
  stopifnot(nrow(proj) >= 1)
  proj <- proj[1, ]
  
  message("Creating RSE (this may download data)...")
  rse <- create_rse(proj)
  
  message("Transforming counts to read counts...")
  assay(rse, "counts") <- transform_counts(rse)
  
  cd <- as.data.frame(colData(rse))
  st_col <- find_col(colnames(cd), c("sample_type"))
  if (is.na(st_col)) stop("Could not find sample_type column in colData.")
  
  cd <- cd %>%
    mutate(sample_type = .data[[st_col]]) %>%
    filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal")) %>%
    mutate(condition = if_else(sample_type == "Primary Tumor", "tumor", "normal"))
  
  rse_sub <- rse[, rownames(cd)]
  
  counts <- assay(rse_sub, "counts")
  counts <- round(counts)
  storage.mode(counts) <- "integer"
  
  coldata <- cd %>%
    transmute(
      sample_id = rownames(cd),
      condition = factor(condition, levels = c("normal", "tumor")),
      sample_type = sample_type
    ) %>%
    as.data.frame()
  rownames(coldata) <- coldata$sample_id
  
  saveRDS(counts, file.path("data/processed", paste0(out_prefix, "_counts.rds")))
  saveRDS(coldata, file.path("data/processed", paste0(out_prefix, "_coldata.rds")))
  
  write_csv(
    as.data.frame(counts) %>% rownames_to_column("gene_id"),
    file.path("data/processed", paste0(out_prefix, "_counts.csv"))
  )
  write_csv(
    coldata %>% rownames_to_column("sample_id"),
    file.path("data/processed", paste0(out_prefix, "_coldata.csv"))
  )
  
  message("Saved: ", out_prefix, " (samples: ", nrow(coldata), ")")
  message("Tumor/Normal: ",
          sum(coldata$condition == "tumor"), "/",
          sum(coldata$condition == "normal"))
  
  invisible(list(counts = counts, coldata = coldata))
}

brca <- prepare_tcga_project("TCGA-BRCA", "tcga_brca")
luad <- prepare_tcga_project("TCGA-LUAD", "tcga_luad")

message("\nAll done.")
