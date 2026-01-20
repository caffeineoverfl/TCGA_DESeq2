# ============================================================
# 01_download_prepare_recount3.R
# Robust TCGA download + preparation via recount3 (schema-flexible)
# Supports TCGA codes like "BRCA" or IDs like "TCGA-BRCA"
# ============================================================

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

cat("=== Script started:", as.character(Sys.time()), "===\n")
cat("Working directory:", getwd(), "\n\n")

normalize_tcga_code <- function(x) {
  x <- toupper(trimws(x))
  sub("^TCGA[-_]", "", x)
}

find_file_recursive <- function(root_dirs, pattern_regex) {
  out <- character(0)
  for (root in root_dirs) {
    root_expanded <- path.expand(root)
    if (!dir.exists(root_expanded)) next
    files <- list.files(
      root_expanded,
      pattern = pattern_regex,
      recursive = TRUE,
      full.names = TRUE,
      include.dirs = FALSE
    )
    if (length(files) > 0) out <- c(out, files)
  }
  unique(out)
}

get_cache_roots <- function() {
  roots <- character(0)
  
  xdg <- Sys.getenv("XDG_CACHE_HOME", unset = NA_character_)
  if (!is.na(xdg) && nzchar(xdg)) roots <- c(roots, xdg)
  
  roots <- c(roots, "~/.cache", "~/.local/share")
  
  rdir <- tryCatch(tools::R_user_dir("recount3", "cache"), error = function(e) NA_character_)
  if (!is.na(rdir) && nzchar(rdir)) roots <- c(roots, rdir)
  
  unique(roots)
}

load_tcga_projects <- function() {
  cat("Fetching recount3 project metadata (this will cache files)...\n")
  invisible(available_projects())
  
  roots <- get_cache_roots()
  cat("Searching for cached TCGA metadata under:\n  - ",
      paste(roots, collapse = "\n  - "), "\n", sep = "")
  
  files <- find_file_recursive(roots, "tcga\\.recount_project\\.MD(\\.gz)?$")
  if (length(files) == 0) {
    stop("Could not find cached file tcga.recount_project.MD(.gz).")
  }
  
  files_gz <- files[str_ends(files, ".gz")]
  candidates <- if (length(files_gz) > 0) files_gz else files
  tcga_path <- candidates[which.max(file.info(candidates)$mtime)]
  
  cat("Found TCGA metadata file:\n  ", tcga_path, "\n", sep = "")
  tcga_projects <- read_tsv(tcga_path, show_col_types = FALSE)
  
  if (!("project" %in% colnames(tcga_projects))) {
    stop("TCGA metadata loaded, but missing required column: project\nColumns: ",
         paste(colnames(tcga_projects), collapse = ", "))
  }
  
  cat("Example projects: ",
      paste(head(sort(unique(tcga_projects$project)), 12), collapse = ", "),
      "\n\n", sep = "")
  
  tcga_projects
}

prepare_tcga_project <- function(tcga_code_or_id, out_prefix, tcga_projects) {
  tcga_code <- normalize_tcga_code(tcga_code_or_id)
  cat("\n=== Preparing TCGA ", tcga_code, " ===\n", sep = "")
  
  proj <- tcga_projects %>% filter(toupper(project) == tcga_code)
  if (nrow(proj) < 1) {
    stop("TCGA code '", tcga_code, "' not found. Available examples:\n",
         paste(head(sort(unique(tcga_projects$project)), 30), collapse = "\n"))
  }
  
  proj1 <- as.data.frame(proj[1, , drop = FALSE])
  
  # Required by create_rse() in your recount3 build
  proj1$project      <- as.character(proj1$project)
  proj1$organism     <- "human"
  proj1$project_home <- "data_sources/tcga"
  if (!("project_type" %in% colnames(proj1))) proj1$project_type <- "data_sources"
  
  cat("Creating RSE (first run may download files)...\n")
  rse <- recount3::create_rse(proj1)
  
  cat("Transforming coverage counts -> read counts...\n")
  assay(rse, "counts") <- recount3::transform_counts(rse)
  
  cd <- as.data.frame(colData(rse))
  st_col <- colnames(cd)[grepl("sample_type", colnames(cd), ignore.case = TRUE)][1]
  if (is.na(st_col)) stop("Could not find 'sample_type' in colData.")
  
  cd <- cd %>%
    mutate(sample_type = .data[[st_col]]) %>%
    filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal")) %>%
    mutate(condition = if_else(sample_type == "Primary Tumor", "tumor", "normal"))
  
  if (nrow(cd) < 2) stop("Too few samples after tumor/normal filtering for ", tcga_code)
  
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
  
  # FIX: coldata already has sample_id column; don't create a duplicate
  write_csv(
    coldata,
    file.path("data/processed", paste0(out_prefix, "_coldata.csv"))
  )
  
  cat("Saved:\n",
      "  data/processed/", out_prefix, "_counts.rds\n",
      "  data/processed/", out_prefix, "_coldata.rds\n", sep = "")
  cat("Samples (tumor/normal): ",
      sum(coldata$condition == "tumor"), "/",
      sum(coldata$condition == "normal"), "\n", sep = "")
  
  invisible(list(counts = counts, coldata = coldata))
}

# ----------------------------
# Main
# ----------------------------

tcga_projects <- load_tcga_projects()

brca <- prepare_tcga_project("BRCA", "tcga_brca", tcga_projects)
luad <- prepare_tcga_project("LUAD", "tcga_luad", tcga_projects)

cat("\n=== Script finished:", as.character(Sys.time()), "===\n")
