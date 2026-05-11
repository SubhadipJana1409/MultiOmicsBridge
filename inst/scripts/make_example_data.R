# ============================================================================
# inst/scripts/make_example_data.R
#
# Creates and saves inst/extdata/example_multiomics.rds — a small,
# self-contained paired multi-omics dataset for use in examples,
# tests, and the package vignette. This script documents the complete
# provenance of the example data, as required by Bioconductor guidelines.
#
# DATA DESCRIPTION:
#   A simulated IBD-inspired paired multi-omics dataset:
#   - 30 samples (15 UC, 15 Control)
#   - 800 host genes (with known IBD-associated signal injected)
#   - 60 microbial taxa (with IBD-associated dysbiosis injected)
#   - Matched sample IDs across both layers
#
# SIGNALS INJECTED (for validation):
#   Host:
#     - Genes 1-15 (CXCL8, S100A8, etc.) upregulated in UC (FC ≈ 4x)
#     - Genes 16-25 (MUC2, TFF3, etc.) downregulated in UC (FC ≈ 0.25x)
#   Microbiome:
#     - Taxa 1-8 (Ruminococcus_gnavus type) enriched in UC (FC ≈ 3x)
#     - Taxa 9-15 (Faecalibacterium_prausnitzii type) depleted in UC (FC ≈ 0.3x)
#
# The cross-omics correlation between CXCL8 (host) and Ruminococcus_gnavus
# (microbiome) is deliberately introduced to test the biomarker network.
#
# OUTPUT:
#   inst/extdata/example_multiomics.rds  (~300 KB, xz-compressed)
#     A named list:
#       $host_se   : SummarizedExperiment (host)
#       $mb_se     : SummarizedExperiment (microbiome)
#       $mae       : MultiAssayExperiment (matched)
#       $outcome   : character vector of condition labels
#       $metadata  : data.frame of sample metadata
#
# TO REGENERATE:
#   Rscript inst/scripts/make_example_data.R
#   # or within an R session:
#   source("inst/scripts/make_example_data.R")
#
# TO LOAD IN EXAMPLES:
#   ex <- readRDS(system.file("extdata","example_multiomics.rds",
#                              package = "MultiOmicsBridge"))
#   mae     <- ex$mae
#   outcome <- ex$outcome
# ============================================================================

cat("================================================================\n")
cat("  MultiOmicsBridge — Build inst/extdata Example Dataset\n")
cat("================================================================\n\n")

suppressPackageStartupMessages({
    library(MultiOmicsBridge)
    library(SummarizedExperiment)
    library(MultiAssayExperiment)
    library(S4Vectors)
})

set.seed(2026L)

# ── Parameters ────────────────────────────────────────────────────────────────

n_samples   <- 30L    # 15 UC + 15 Control
n_genes     <- 800L
n_taxa      <- 60L
n_uc        <- 15L
n_ctrl      <- 15L

# Injected signal indices
up_gene_idx   <- 1:15     # upregulated in UC
down_gene_idx <- 16:25    # downregulated in UC
up_taxa_idx   <- 1:8      # enriched in UC
down_taxa_idx <- 9:15     # depleted in UC

sample_ids <- paste0("Sample", seq_len(n_samples))
outcome    <- rep(c("UC", "Control"), each = c(n_uc, n_ctrl))
is_uc      <- outcome == "UC"

cat(sprintf("Generating example dataset:\n"))
cat(sprintf("  Samples : %d (%d UC, %d Control)\n", n_samples, n_uc, n_ctrl))
cat(sprintf("  Genes   : %d\n", n_genes))
cat(sprintf("  Taxa    : %d\n", n_taxa))

# ── Host gene names (real-sounding IBD-relevant names) ────────────────────────

# Named signal genes
ibd_up_names   <- c("CXCL8","S100A8","S100A9","IL6","MMP9",
                     "STAT3","OSM","FCGR3A","ITGAM","CXCL10",
                     "IL1B","TNF","S100A12","LCN2","CXCL1")
ibd_down_names <- c("MUC2","TFF3","PPARG","FABP1","FABP4",
                     "HMGCS2","SLC5A1","CLDN8","ATOH1","SPDEF")

# Fill rest with neutral names
n_neutral_genes <- n_genes - length(ibd_up_names) - length(ibd_down_names)
neutral_genes   <- paste0("GENE_", sprintf("%04d", seq_len(n_neutral_genes)))
all_gene_names  <- c(ibd_up_names, ibd_down_names, neutral_genes)

# ── Microbial taxon names ─────────────────────────────────────────────────────

uc_up_taxa_names <- c("Ruminococcus_gnavus","Veillonella_parvula",
                       "Haemophilus_parainfluenzae","Escherichia_coli",
                       "Streptococcus_salivarius","Clostridium_ramosum",
                       "Fusobacterium_nucleatum","Klebsiella_pneumoniae")
uc_dn_taxa_names <- c("Faecalibacterium_prausnitzii","Roseburia_intestinalis",
                       "Blautia_obeum","Coprococcus_comes",
                       "Akkermansia_muciniphila","Bifidobacterium_longum",
                       "Eubacterium_hallii")
n_neutral_taxa   <- n_taxa - length(uc_up_taxa_names) - length(uc_dn_taxa_names)
neutral_taxa     <- paste0("uncultured_bacterium_sp", sprintf("%03d",
                            seq_len(n_neutral_taxa)))
all_taxa_names   <- c(uc_up_taxa_names, uc_dn_taxa_names, neutral_taxa)

# ── Generate host count matrix ────────────────────────────────────────────────

cat("\nGenerating host RNA-seq counts...\n")

host_counts <- matrix(
    rnbinom(n_genes * n_samples, mu = 180, size = 3),
    nrow = n_genes, ncol = n_samples,
    dimnames = list(all_gene_names, sample_ids)
)

# Apply UC-specific fold changes with sample-level noise
for (i in up_gene_idx) {
    fc_uc <- rlnorm(n_uc,   meanlog = log(4.0), sdlog = 0.35)
    host_counts[i, is_uc]  <- round(host_counts[i, is_uc] * fc_uc)
}
for (i in down_gene_idx) {
    fc_uc <- rlnorm(n_uc, meanlog = log(0.25), sdlog = 0.30)
    host_counts[i, is_uc] <- pmax(1L, round(host_counts[i, is_uc] * fc_uc))
}

# Simulate realistic library size variation (50-200M reads range)
lib_scale <- rlnorm(n_samples, meanlog = 0, sdlog = 0.3)
for (j in seq_len(n_samples)) {
    host_counts[, j] <- round(host_counts[, j] * lib_scale[j])
}
storage.mode(host_counts) <- "integer"

# ── Generate microbiome count matrix ─────────────────────────────────────────

cat("Generating microbiome counts...\n")

mb_counts <- matrix(
    rnbinom(n_taxa * n_samples, mu = 600, size = 2),
    nrow = n_taxa, ncol = n_samples,
    dimnames = list(all_taxa_names, sample_ids)
)

for (i in up_taxa_idx) {
    fc_uc <- rlnorm(n_uc, meanlog = log(3.5), sdlog = 0.4)
    mb_counts[i, is_uc] <- round(mb_counts[i, is_uc] * fc_uc)
}
for (i in down_taxa_idx) {
    fc_uc <- rlnorm(n_uc, meanlog = log(0.3), sdlog = 0.35)
    mb_counts[i, is_uc] <- pmax(1L, round(mb_counts[i, is_uc] * fc_uc))
}

# Inject deliberate CXCL8 ↔ Ruminococcus_gnavus cross-omics correlation
# This should appear in the biomarker network heatmap
rg_profile <- mb_counts["Ruminococcus_gnavus", ] / max(mb_counts["Ruminococcus_gnavus", ])
cxcl8_base  <- host_counts["CXCL8", ]
host_counts["CXCL8", ] <- round(cxcl8_base * 0.4 + cxcl8_base * 0.6 *
                                   rg_profile * max(cxcl8_base) /
                                   max(rg_profile * max(cxcl8_base)))

storage.mode(mb_counts) <- "integer"

# ── Sample metadata ───────────────────────────────────────────────────────────

metadata <- data.frame(
    condition  = outcome,
    age        = sample(25:70, n_samples, replace = TRUE),
    sex        = sample(c("M", "F"), n_samples, replace = TRUE),
    bmi        = round(rnorm(n_samples, mean = 24, sd = 3.5), 1),
    disease_duration_yrs = ifelse(is_uc,
                                   round(rlnorm(n_samples, 1.5, 0.6), 1), NA),
    row.names  = sample_ids,
    stringsAsFactors = FALSE
)

# ── Build SummarizedExperiments and MAE ──────────────────────────────────────

cat("Building SummarizedExperiment objects...\n")

host_se <- loadHostData(
    counts   = host_counts,
    col_data = metadata,
    min_count = 1L
)

mb_se <- loadMicrobiomeData(
    taxa_table     = mb_counts,
    col_data       = metadata,
    normalization  = "CLR",
    pseudocount    = 0.5,
    min_prevalence = 0
)

mae <- matchSamples(host_se, mb_se, min_paired = 10)

# ── Quick validation ──────────────────────────────────────────────────────────

cat("\nValidating example data...\n")

# Verify signal recovery with a quick pipeline run
dr_test <- jointDimReduction(mae, outcome,
                              n_components = 2,
                              n_features_host = 20,
                              n_features_mb = 10)

bm_test <- biomarkerDiscovery(mae, dr_test, n_biomarkers = 15)
bm_df   <- as.data.frame(bm_test)

host_bm_test <- bm_df$feature[bm_df$omics_layer == "host"]
mb_bm_test   <- bm_df$feature[bm_df$omics_layer == "microbiome"]

up_rec   <- intersect(host_bm_test, ibd_up_names)
down_rec <- intersect(host_bm_test, ibd_down_names)
taxa_rec <- intersect(mb_bm_test, uc_up_taxa_names)

cat(sprintf("  IBD-up genes recovered   : %d / %d\n",
            length(up_rec), length(ibd_up_names)))
cat(sprintf("  IBD-down genes recovered : %d / %d\n",
            length(down_rec), length(ibd_down_names)))
cat(sprintf("  IBD taxa recovered       : %d / %d\n",
            length(taxa_rec), length(uc_up_taxa_names)))

# CXCL8–Ruminococcus correlation
cxcl8_row    <- bm_df[bm_df$feature == "CXCL8", ]
rg_in_bm     <- "Ruminococcus_gnavus" %in% mb_bm_test
cxcl8_in_bm  <- "CXCL8" %in% host_bm_test

if (!is.null(cxcl8_row$top_partner) && nrow(cxcl8_row) > 0L) {
    cat(sprintf("  CXCL8 top partner        : %s (r=%.3f)\n",
                cxcl8_row$top_partner[1],
                cxcl8_row$max_cross_cor[1]))
}
cat(sprintf("  CXCL8 in biomarkers      : %s\n", cxcl8_in_bm))
cat(sprintf("  Ruminococcus in biomarkers: %s\n", rg_in_bm))

if (length(up_rec) == 0L)
    warning("No IBD-up genes recovered. Check signal injection.")

# ── Save example data ─────────────────────────────────────────────────────────

outdir   <- file.path("inst", "extdata")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(outdir, "example_multiomics.rds")

example_data <- list(
    host_se  = host_se,
    mb_se    = mb_se,
    mae      = mae,
    outcome  = outcome,
    metadata = metadata,
    # Ground truth for test validation
    ground_truth = list(
        ibd_up_genes   = ibd_up_names,
        ibd_down_genes = ibd_down_names,
        uc_up_taxa     = uc_up_taxa_names,
        uc_down_taxa   = uc_dn_taxa_names,
        correlated_pair = c(host = "CXCL8",
                             microbiome = "Ruminococcus_gnavus")
    )
)

saveRDS(example_data, file = out_path, compress = "xz")

fsize <- file.size(out_path)
cat(sprintf("\n=== Example dataset saved ===\n"))
cat(sprintf("  Path       : %s\n", out_path))
cat(sprintf("  Size       : %.1f KB (xz-compressed)\n", fsize / 1024))
cat(sprintf("  Samples    : %d (%d UC, %d Control)\n", n_samples, n_uc, n_ctrl))
cat(sprintf("  Host genes : %d\n", nrow(host_se)))
cat(sprintf("  Taxa       : %d\n", nrow(mb_se)))
cat(sprintf("  Assays     : host=%s | mb=%s\n",
            paste(assayNames(host_se), collapse = "+"),
            paste(assayNames(mb_se),   collapse = "+")))

cat("\nTo load in your R session:\n")
cat('  ex <- readRDS(system.file("extdata", "example_multiomics.rds",\n')
cat('                              package = "MultiOmicsBridge"))\n')
cat('  mae     <- ex$mae\n')
cat('  outcome <- ex$outcome\n')
cat('\nGeneration complete.\n')

sessionInfo()
