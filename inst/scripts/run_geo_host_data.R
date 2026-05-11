# ============================================================================
# inst/scripts/run_geo_host_data.R
#
# How to substitute REAL GEO host RNA-seq data into MultiOmicsBridge.
# Uses GEOquery to download GSE87466 (UC vs Control rectal RNA-seq,
# 90 samples, Vanhove et al. 2018) as the host transcriptomics layer.
#
# DATASET:
#   GSE87466 — RNA-seq from rectal biopsies (45 UC, 45 control)
#   Vanhove W et al. (2018) Gut 67:1573-1581
#
# MICROBIOME:
#   FranzsosaEA_2019 from curatedMetagenomicData (same as run_ibdmdb_demo.R)
#   NOTE: Not truly paired (different individuals), but condition-matched.
#         This demonstrates the workflow for condition-matched multi-omics.
#         For truly paired samples, use the full IBDMDB portal download
#         (see run_ibdmdb_full_paired.R).
#
# REQUIRES:
#   BiocManager::install(c("GEOquery", "curatedMetagenomicData"))
#
# RUNTIME: ~15-25 minutes (GEO download + processing)
#
# USAGE:
#   source("inst/scripts/run_geo_host_data.R")
# ============================================================================

cat("==========================================================\n")
cat("  MultiOmicsBridge — Real GEO Host RNA-seq Demo\n")
cat("  GSE87466 (UC RNA-seq) + FranzsosaEA_2019 (Microbiome)\n")
cat("==========================================================\n\n")

# ── 0. Package checks ─────────────────────────────────────────────────────────

required_pkgs <- c(
    "MultiOmicsBridge",
    "GEOquery",
    "curatedMetagenomicData",
    "SummarizedExperiment",
    "MultiAssayExperiment",
    "limma", "edgeR"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs,
    requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(missing_pkgs, ask = FALSE)
}

suppressPackageStartupMessages({
    library(MultiOmicsBridge)
    library(GEOquery)
    library(curatedMetagenomicData)
    library(SummarizedExperiment)
    library(MultiAssayExperiment)
})

outdir <- file.path(tempdir(), "MultiOmicsBridge_GEO")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Output directory: %s\n\n", outdir))

t0 <- Sys.time()

# ── 1. Download GSE87466 from GEO ─────────────────────────────────────────────

cat("=== STEP 1: Download GSE87466 from GEO ===\n")
cat("Dataset: Vanhove W et al. (2018) Gut | UC vs Control RNA-seq\n")
cat("(Downloading from GEO — may take a few minutes)\n\n")

# Download expression data
gse <- tryCatch(
    GEOquery::getGEO("GSE87466",
                     GSEMatrix    = TRUE,
                     getGPL       = FALSE,
                     destdir      = outdir),
    error = function(e) {
        stop("GEO download failed: ", conditionMessage(e),
             "\nCheck your internet connection and try again.")
    }
)

# getGEO returns a list; take the first element
gse_obj <- gse[[1]]
cat(sprintf("Downloaded %d samples with %d features\n",
            ncol(exprs(gse_obj)), nrow(exprs(gse_obj))))

# ── 2. Extract and process expression data ────────────────────────────────────

cat("\n=== STEP 2: Process GEO Expression Data ===\n")

# Extract phenotype data
pdata <- pData(gse_obj)
cat("Sample metadata columns:\n")
print(names(pdata))

# Identify condition column (look for "characteristics_ch1" columns)
char_cols <- grep("characteristics_ch1", names(pdata), value = TRUE)
cat("\nCharacteristic columns:\n")
print(pdata[1:3, char_cols, drop = FALSE])

# GSE87466 specific: condition is in characteristics_ch1
# "diagnosis: UC" or "diagnosis: normal"
cond_col_geo <- char_cols[1]
conditions_raw <- as.character(pdata[[cond_col_geo]])

# Parse condition labels
# Parse condition labels using base R (no dplyr dependency)
condition_geo <- ifelse(
    grepl("UC|ulcerative", conditions_raw, ignore.case = TRUE), "UC",
    ifelse(grepl("normal|control|healthy", conditions_raw,
                 ignore.case = TRUE), "Control", NA_character_)
)

cat("\nCondition distribution:\n")
print(table(condition_geo, useNA = "always"))

# Keep only UC and Control
keep_geo <- !is.na(condition_geo)
gse_sub  <- gse_obj[, keep_geo]
cond_sub <- condition_geo[keep_geo]

cat(sprintf("\nRetaining: %d samples (%d UC, %d Control)\n",
            sum(keep_geo), sum(cond_sub == "UC"), sum(cond_sub == "Control")))

# Extract expression matrix
# GSE87466 is a microarray dataset (Affymetrix) — convert to count-like
# matrix by taking 2^x to undo log2 transformation and rounding
expr_mat <- exprs(gse_sub)

# Check if data is log-transformed (typical for GEO processed data)
if (max(expr_mat, na.rm = TRUE) < 30) {
    cat("Data appears log2-transformed. Converting to linear scale.\n")
    expr_mat <- 2^expr_mat
}

# For GSE87466 specifically, if it's an RNA-seq dataset with count data:
# Check the range to determine data type
data_range <- range(expr_mat, na.rm = TRUE)
cat(sprintf("Expression range: %.1f to %.1f\n", data_range[1], data_range[2]))

# Round to integer counts for loadHostData (safe for both microarray and RNA-seq)
expr_counts <- round(expr_mat)
expr_counts[is.na(expr_counts)] <- 0L
expr_counts[expr_counts < 0]    <- 0L
storage.mode(expr_counts) <- "integer"

# Generate stable sample names
new_sample_names <- paste0("GEO_S", seq_len(ncol(expr_counts)))
colnames(expr_counts) <- new_sample_names

geo_col_data <- data.frame(
    condition  = cond_sub,
    geo_sample = colnames(gse_sub),
    row.names  = new_sample_names,
    stringsAsFactors = FALSE
)

cat(sprintf("\nLoading %d genes × %d samples into MultiOmicsBridge...\n",
            nrow(expr_counts), ncol(expr_counts)))

host_se_geo <- loadHostData(
    counts    = expr_counts,
    col_data  = geo_col_data,
    min_count = 5L
)
cat("Host SummarizedExperiment (REAL GEO DATA):\n")
print(host_se_geo)

# ── 3. Download FranzsosaEA_2019 microbiome data ──────────────────────────────

cat("\n=== STEP 3: Download IBD Microbiome (curatedMetagenomicData) ===\n")

tse_list <- curatedMetagenomicData(
    "FranzsosaEA_2019.relative_abundance",
    dryrun  = FALSE,
    rownames = "short"
)
tse <- tse_list[[1]]
sample_meta <- as.data.frame(colData(tse))

# Filter to UC vs nonIBD
cond_col_mb <- intersect(c("disease","study_condition","condition"),
                          names(sample_meta))[1]
keep_mb <- sample_meta[[cond_col_mb]] %in% c("UC", "nonIBD")
tse_mb  <- tse[, keep_mb]
meta_mb <- sample_meta[keep_mb, ]

mb_mat_raw      <- as.matrix(assay(tse_mb, "relative_abundance"))
mb_pseudocounts <- round(mb_mat_raw * 1e6)
storage.mode(mb_pseudocounts) <- "integer"

# Prevalence filter
prev <- rowMeans(mb_pseudocounts > 0)
mb_pseudocounts <- mb_pseudocounts[prev >= 0.20, ]

# Rename samples for matching (condition-matched, not individually paired)
n_mb_uc   <- sum(meta_mb[[cond_col_mb]] == "UC")
n_mb_ctrl <- sum(meta_mb[[cond_col_mb]] == "nonIBD")

# ── 4. Balance datasets for condition-matched integration ─────────────────────

cat("\n=== STEP 4: Balance Datasets for Condition-Matched Integration ===\n")
cat("NOTE: These are condition-matched (not individual-paired) samples.\n")
cat("      UC patients from GSE87466 are matched to UC from IBDMDB,\n")
cat("      controls from GSE87466 are matched to nonIBD from IBDMDB.\n\n")

# Find the minimum number of samples in each class across both datasets
n_geo_uc   <- sum(cond_sub == "UC")
n_geo_ctrl <- sum(cond_sub == "Control")

n_use_uc   <- min(n_geo_uc, n_mb_uc)
n_use_ctrl <- min(n_geo_ctrl, n_mb_ctrl)

cat(sprintf("UC samples   — GEO: %d, Microbiome: %d → using %d\n",
            n_geo_uc, n_mb_uc, n_use_uc))
cat(sprintf("Control samples — GEO: %d, Microbiome: %d → using %d\n",
            n_geo_ctrl, n_mb_ctrl, n_use_ctrl))

# Select balanced subsets
set.seed(42)
geo_uc_idx   <- which(cond_sub == "UC")
geo_ctrl_idx <- which(cond_sub == "Control")
mb_uc_idx    <- which(meta_mb[[cond_col_mb]] == "UC")
mb_ctrl_idx  <- which(meta_mb[[cond_col_mb]] == "nonIBD")

geo_sel <- c(
    sample(geo_uc_idx,   n_use_uc),
    sample(geo_ctrl_idx, n_use_ctrl)
)
mb_sel <- c(
    sample(mb_uc_idx,   n_use_uc),
    sample(mb_ctrl_idx, n_use_ctrl)
)

n_total <- n_use_uc + n_use_ctrl
shared_names <- paste0("Sample", seq_len(n_total))
shared_cond  <- c(rep("UC", n_use_uc), rep("Control", n_use_ctrl))

# Subset and rename
host_mat_sel <- assay(host_se_geo, "counts")[, geo_sel, drop = FALSE]
colnames(host_mat_sel) <- shared_names
host_cd <- data.frame(condition = shared_cond, row.names = shared_names)

mb_mat_sel <- mb_pseudocounts[, mb_sel, drop = FALSE]
colnames(mb_mat_sel) <- shared_names
mb_cd <- data.frame(condition = shared_cond, row.names = shared_names)

# Reload with matched names
host_se_matched <- loadHostData(host_mat_sel, col_data = host_cd, min_count = 5L)
mb_se_matched   <- loadMicrobiomeData(mb_mat_sel, col_data = mb_cd,
                                       normalization = "CLR",
                                       min_prevalence = 0)

# ── 5. Match samples and run pipeline ────────────────────────────────────────

cat("\n=== STEP 5: Match Samples and Run Full Pipeline ===\n")

mae_geo <- matchSamples(host_se_matched, mb_se_matched, min_paired = 10)
outcome_geo <- colData(mae_geo)$condition

result_geo <- MultiOmicsBridgeAnalysis(
    mae_geo,
    outcome         = outcome_geo,
    n_components    = 2L,
    n_features_host = 50L,
    n_features_mb   = 20L,
    n_biomarkers    = 40L,
    cv_folds        = 5L,
    seed            = 42L
)

# ── 6. Results ────────────────────────────────────────────────────────────────

cat("\n=== STEP 6: Results ===\n")
show(result_geo)

cr_geo <- performance(result_geo)
cat("\nClassifier AUC (5-fold CV):\n")
cat(sprintf("  Host only      : %.3f ± %.3f\n",
            cr_geo$host_only$mean_auc, cr_geo$host_only$sd_auc))
cat(sprintf("  Microbiome only: %.3f ± %.3f\n",
            cr_geo$microbiome_only$mean_auc, cr_geo$microbiome_only$sd_auc))
cat(sprintf("  Joint          : %.3f ± %.3f\n",
            cr_geo$joint$mean_auc, cr_geo$joint$sd_auc))

# Plots
cat("\nGenerating plots...\n")
p1 <- plotIntegration(result_geo, outcome = outcome_geo,
                       colours = c(UC = "#E24B4A", Control = "#378ADD"))
ggsave(file.path(outdir, "geo_integration.png"), p1, width = 9, height = 6)

p2 <- plotClassifierComparison(result_geo, type = "bar")
ggsave(file.path(outdir, "geo_classifier.png"), p2, width = 7, height = 5)

p3 <- plotBiomarkerNetwork(result_geo, mae_geo, n_host = 15, n_mb = 10)
ggsave(file.path(outdir, "geo_network.png"), p3, width = 10, height = 8)

cat(sprintf("\nAll plots saved to: %s\n", outdir))

generateReport(result_geo, n_top = 10)

elapsed_total <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("\nTotal runtime: %.1f seconds\n", elapsed_total))
cat("\n=== GEO Real Data Demo Complete ===\n")
sessionInfo()
