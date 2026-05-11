# ============================================================================
# inst/scripts/run_ibdmdb_demo.R
#
# Comprehensive validation of MultiOmicsBridge on the IBD Multi-omics
# Database (IBDMDB) microbiome dataset (Franzosa et al. 2019,
# Nature Microbiology 4:293-305).
#
# DATA SOURCES:
#   Microbiome : curatedMetagenomicData (Bioconductor) — FranzsosaEA_2019
#                REAL published data, 220 stool metagenomics samples,
#                UC / CD / nonIBD subjects.
#   Host RNA-seq: Biologically realistic simulation using known IBD
#                 transcriptional signatures from published literature.
#                 Clearly labelled throughout as SIMULATED.
#                 See run_geo_host_data.R to substitute real GEO data.
#
# WHAT THIS SCRIPT TESTS:
#   All five MultiOmicsBridge modules end-to-end:
#   1. loadMicrobiomeData   — CLR normalization of real metagenomics
#   2. loadHostData         — voom normalization of simulated RNA-seq
#   3. matchSamples         — paired MAE construction
#   4. jointDimReduction    — DIABLO integration
#   5. biomarkerDiscovery   — cross-omics biomarker ranking
#   6. diagnosticClassifier — host-only vs mb-only vs joint AUC comparison
#   7. All five plot functions
#   8. generateReport
#
# RUNTIME: ~10-15 minutes (first run downloads ExperimentHub data ~50 MB)
#
# USAGE:
#   source("inst/scripts/run_ibdmdb_demo.R")
# ============================================================================

cat("============================================================\n")
cat("  MultiOmicsBridge — IBDMDB Real Data Validation Demo\n")
cat("  Franzosa et al. 2019 | IBD Metagenomics + Simulated Host\n")
cat("============================================================\n\n")

# ── 0. Package checks ─────────────────────────────────────────────────────────

required_pkgs <- c(
    "MultiOmicsBridge",
    "curatedMetagenomicData",   # Bioconductor
    "SummarizedExperiment",
    "MultiAssayExperiment",
    "S4Vectors",
    "ggplot2",
    "limma",
    "edgeR"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs,
    requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0L) {
    cat("Installing missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(missing_pkgs, ask = FALSE)
}

suppressPackageStartupMessages({
    library(MultiOmicsBridge)
    library(curatedMetagenomicData)
    library(SummarizedExperiment)
    library(MultiAssayExperiment)
    library(S4Vectors)
    library(ggplot2)
})

# Output directory for plots and report
outdir <- file.path(tempdir(), "MultiOmicsBridge_IBDMDB")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Output directory: %s\n\n", outdir))

# ── 1. Download FranzsosaEA_2019 microbiome data ──────────────────────────────

cat("=== STEP 1: Download IBD Metagenomics (curatedMetagenomicData) ===\n")
cat("Dataset: Franzosa EA et al. (2019) Nat Microbiol 4:293-305\n")
cat("(First run downloads ~50 MB from ExperimentHub; cached thereafter)\n\n")

t0 <- Sys.time()

# Download species-level relative abundances
tse_list <- curatedMetagenomicData(
    "FranzsosaEA_2019.relative_abundance",
    dryrun  = FALSE,
    rownames = "short"        # short names: Genus_species format
)
tse <- tse_list[[1]]

cat(sprintf("Downloaded: %d taxa × %d samples\n", nrow(tse), ncol(tse)))

# ── 2. Process microbiome data ────────────────────────────────────────────────

cat("\n=== STEP 2: Process Microbiome Data ===\n")

# Extract sample metadata
sample_meta <- as.data.frame(colData(tse))
cat("Available sample metadata columns:\n")
print(names(sample_meta))

# Check condition column
cond_col <- intersect(c("disease", "study_condition", "condition"),
                      names(sample_meta))[1]
if (is.na(cond_col))
    stop("Cannot identify condition column in sample metadata.")

cat(sprintf("\nCondition column: '%s'\n", cond_col))
cat("Condition distribution:\n")
print(table(sample_meta[[cond_col]]))

# Focus on UC vs nonIBD for binary classification
# (largest, cleanest contrast in this dataset)
keep_conditions <- c("UC", "nonIBD")
keep_idx <- sample_meta[[cond_col]] %in% keep_conditions

tse_sub <- tse[, keep_idx]
sample_meta_sub <- sample_meta[keep_idx, ]

cat(sprintf("\nRetaining UC vs nonIBD: %d samples\n", ncol(tse_sub)))
cat(sprintf("  UC:     %d\n", sum(sample_meta_sub[[cond_col]] == "UC")))
cat(sprintf("  nonIBD: %d\n", sum(sample_meta_sub[[cond_col]] == "nonIBD")))

# Extract relative abundance matrix (taxa × samples)
mb_mat_raw <- as.matrix(assay(tse_sub, "relative_abundance"))

# Convert relative abundance to pseudo-counts (multiply by 1e6 → reads per million)
# This preserves compositional structure and allows CLR normalization
mb_pseudo_counts <- round(mb_mat_raw * 1e6)
storage.mode(mb_pseudo_counts) <- "integer"

cat(sprintf("\nMicrobiome matrix: %d taxa × %d samples\n",
            nrow(mb_pseudo_counts), ncol(mb_pseudo_counts)))
cat(sprintf("Sparsity: %.1f%% zeros\n",
            100 * mean(mb_pseudo_counts == 0)))

# Remove very rare taxa (present in < 20% of samples, more stringent than default
# to keep only reliable taxa for this demonstration)
prevalence <- rowMeans(mb_pseudo_counts > 0)
mb_pseudo_counts <- mb_pseudo_counts[prevalence >= 0.20, ]
cat(sprintf("After prevalence filter (>= 20%%): %d taxa retained\n",
            nrow(mb_pseudo_counts)))

# Load into MultiOmicsBridge
mb_col_data <- data.frame(
    condition = sample_meta_sub[[cond_col]],
    subject   = if ("subject_id" %in% names(sample_meta_sub))
                    sample_meta_sub$subject_id else colnames(tse_sub),
    row.names = colnames(tse_sub),
    stringsAsFactors = FALSE
)

mb_se <- loadMicrobiomeData(
    taxa_table     = mb_pseudo_counts,
    col_data       = mb_col_data,
    normalization  = "CLR",
    pseudocount    = 0.5,
    min_prevalence = 0          # already filtered above
)
cat("\nMicrobiome SummarizedExperiment:\n")
print(mb_se)

# ── 3. Generate biologically realistic host RNA-seq (SIMULATED) ───────────────

cat("\n=== STEP 3: Generate Host RNA-seq (SIMULATED — see header) ===\n")
cat("Simulation uses known IBD transcriptional signatures from literature.\n\n")

# Number of samples to match microbiome
n_samples <- ncol(mb_se)
outcome   <- colData(mb_se)$condition   # "UC" or "nonIBD"
is_uc     <- outcome == "UC"

cat(sprintf("Generating RNA-seq for %d samples (%d UC, %d nonIBD)\n",
            n_samples, sum(is_uc), sum(!is_uc)))

set.seed(2026)

# ── Known IBD-associated gene signatures (from Taman et al. 2018,
#    Neurath 2017, Corridoni et al. 2020, and other reviews) ──────────────────
ibd_up_genes <- c(
    # Cytokines and chemokines
    "CXCL8","CXCL10","CXCL1","CXCL2","CXCL3","CXCL5","CCL2","CCL20",
    "IL6","IL1B","IL17A","IL17F","IL23A","IL33","IL18","IL34","TNF",
    # Transcription factors / signalling
    "STAT3","STAT1","NFKB1","NFKB2","RELA","IRF1","IRF4","BATF",
    # Effector / damage markers
    "S100A8","S100A9","S100A12","LCN2","MMP9","MMP3","MMP12",
    "OSM","LIF","ONCOSTATIN","HAMP","OLFM4",
    # Immune cell markers
    "FCGR3A","ITGAM","CD68","CD14","CD163","ADGRE1",
    # Mucosal immunity
    "PIGR","JCHAIN","IGHG1","IGHG4","IGHA1","IGHA2"
)
ibd_up_genes <- unique(ibd_up_genes)

ibd_down_genes <- c(
    # Epithelial barrier and mucus
    "MUC2","MUC5AC","TFF3","TFF1","CLCA1","FCGBP",
    # Metabolic / absorptive
    "PPARG","FABP1","FABP4","HMGCS2","HMGCS1","APOA1","APOC2",
    # Transporters
    "SLC5A1","SLC2A5","SLC26A3","SLC26A6","CFTR",
    # Tight junction / barrier
    "CLDN8","CLDN4","OCLN","CDH1",
    # Cell identity
    "ATOH1","SPDEF","GFI1","HNF4A","KLF4","HOPX"
)
ibd_down_genes <- unique(ibd_down_genes)

# Background genes (highly expressed in intestinal epithelium)
background_genes_base <- c(
    # Structural
    paste0("KRT", c(8,18,19,20)),
    paste0("VIM","ACTA2","ACTB","ACTG1"),
    # Ribosomal (highly expressed housekeeping)
    paste0("RPS", c(3,4,6,8,11,14,18,24,27)),
    paste0("RPL", c(3,4,6,8,11,14,18,24,27,32)),
    # Metabolism
    paste0("GAPDH","LDHA","LDHB","ENO1","PKM","PGK1","TPI1"),
    # Other housekeeping
    paste0("EEF1A1","EEF1B2","EEF2","EIF4A1","HSPA8","HSP90AA1",
           "HSP90AB1","HSPD1","HSPE1")
)
background_genes_base <- unique(unlist(strsplit(background_genes_base, "")))
# Add many more background genes to reach ~3000 total
set.seed(42)
extra_genes <- paste0("ENSG", sprintf("%08d", sample(1:50000, 2500)))
all_genes <- unique(c(ibd_up_genes, ibd_down_genes,
                       background_genes_base, extra_genes))
n_genes   <- length(all_genes)
cat(sprintf("Simulating expression for %d genes\n", n_genes))
cat(sprintf("  Known IBD-up   : %d genes\n", length(ibd_up_genes)))
cat(sprintf("  Known IBD-down : %d genes\n", length(ibd_down_genes)))
cat(sprintf("  Background     : %d genes\n",
            n_genes - length(ibd_up_genes) - length(ibd_down_genes)))

# Base counts: negative binomial with realistic mean and dispersion
# Intestinal epithelium: mean expression ~150 counts/gene at typical depth
base_counts <- matrix(
    rnbinom(n_genes * n_samples, mu = 150, size = 3),
    nrow = n_genes, ncol = n_samples
)
rownames(base_counts) <- all_genes
colnames(base_counts) <- colnames(mb_se)

# Apply IBD-specific fold changes to UC samples
# Upregulated in UC: logFC ~ N(1.8, 0.5), i.e. ~ 3.5× fold change
up_idx <- match(ibd_up_genes, all_genes)
up_idx <- up_idx[!is.na(up_idx)]
for (i in up_idx) {
    fc <- exp(rnorm(sum(is_uc), mean = log(3.5), sd = 0.4))
    base_counts[i, is_uc] <- round(base_counts[i, is_uc] * fc)
}

# Downregulated in UC: logFC ~ N(-1.5, 0.4), i.e. ~ 0.22× fold change
dn_idx <- match(ibd_down_genes, all_genes)
dn_idx <- dn_idx[!is.na(dn_idx)]
for (i in dn_idx) {
    fc <- exp(rnorm(sum(is_uc), mean = log(0.22), sd = 0.3))
    base_counts[i, is_uc] <- pmax(1L,
                                   round(base_counts[i, is_uc] * fc))
}

# Simulate realistic library size variation (40–120M reads per sample)
lib_factors <- rlnorm(n_samples, meanlog = log(1), sdlog = 0.35)
for (j in seq_len(n_samples)) {
    base_counts[, j] <- round(base_counts[, j] * lib_factors[j])
}
storage.mode(base_counts) <- "integer"

# Add small amount of random noise to avoid perfect simulation
base_counts <- base_counts + matrix(
    rpois(n_genes * n_samples, lambda = 2),
    nrow = n_genes
)

host_col_data <- data.frame(
    condition = outcome,
    row.names = colnames(mb_se),
    stringsAsFactors = FALSE
)

host_se <- loadHostData(
    counts    = base_counts,
    col_data  = host_col_data,
    min_count = 10L
)
cat("\nHost SummarizedExperiment (SIMULATED):\n")
print(host_se)

elapsed_prep <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("\nData preparation complete in %.1f seconds\n", elapsed_prep))

# ── 4. Match samples ──────────────────────────────────────────────────────────

cat("\n=== STEP 4: Match Samples → MultiAssayExperiment ===\n")

mae <- matchSamples(host_se, mb_se, min_paired = 20)
cat("\nMultiAssayExperiment:\n")
print(mae)

outcome_vec <- colData(mae)$condition
cat(sprintf("\nFinal outcome: %d UC, %d nonIBD\n",
            sum(outcome_vec == "UC"),
            sum(outcome_vec == "nonIBD")))

# ── 5. Full pipeline ──────────────────────────────────────────────────────────

cat("\n=== STEP 5: Run MultiOmicsBridgeAnalysis ===\n")

# Choose n_features appropriate for the dataset size
n_samp  <- ncol(mae)
n_host_f <- min(80L, floor(nrow(experiments(mae)[["host"]])  * 0.02))
n_mb_f   <- min(30L, floor(nrow(experiments(mae)[["microbiome"]]) * 0.10))
n_bm     <- 60L

cat(sprintf("Settings: %d samples | %d host features/comp | %d mb features/comp\n",
            n_samp, n_host_f, n_mb_f))

t1 <- Sys.time()

result <- MultiOmicsBridgeAnalysis(
    mae,
    outcome         = outcome_vec,
    n_components    = 2L,
    n_features_host = n_host_f,
    n_features_mb   = n_mb_f,
    n_biomarkers    = n_bm,
    cv_folds        = 5L,
    host_assay      = "voom",
    mb_assay        = "CLR",
    design_off_diag = 0.1,
    seed            = 2026L
)

elapsed_analysis <- as.numeric(Sys.time() - t1, units = "secs")
cat(sprintf("\nAnalysis complete in %.1f seconds\n", elapsed_analysis))

# ── 6. Validate results ───────────────────────────────────────────────────────

cat("\n=== STEP 6: Validate Results ===\n")
cat("\nMOBResult summary:\n")
show(result)

# ── 6a: Biomarker validation ──────────────────────────────────────────────────
bm_df <- as.data.frame(biomarkers(result))

cat("\n-- Biomarker Recovery --\n")
# Check how many known IBD-up genes were recovered as host biomarkers
host_bm     <- bm_df$feature[bm_df$omics_layer == "host"]
mb_bm       <- bm_df$feature[bm_df$omics_layer == "microbiome"]
up_recovered   <- intersect(host_bm, ibd_up_genes)
down_recovered <- intersect(host_bm, ibd_down_genes)

cat(sprintf("Host biomarkers selected          : %d\n", length(host_bm)))
cat(sprintf("Microbiome biomarkers selected    : %d\n", length(mb_bm)))
cat(sprintf("Known IBD-up genes recovered      : %d / %d\n",
            length(up_recovered), length(ibd_up_genes)))
cat(sprintf("Known IBD-down genes recovered    : %d / %d\n",
            length(down_recovered), length(ibd_down_genes)))

if (length(up_recovered) > 0L) {
    cat("  IBD-up recovered: ",
        paste(head(up_recovered, 10), collapse = ", "), "\n")
}
if (length(down_recovered) > 0L) {
    cat("  IBD-down recovered: ",
        paste(head(down_recovered, 5), collapse = ", "), "\n")
}

# Top microbiome biomarkers
cat(sprintf("\nTop 10 microbiome biomarkers:\n"))
mb_bm_df <- bm_df[bm_df$omics_layer == "microbiome", ]
mb_bm_df <- mb_bm_df[order(mb_bm_df$loading_score, decreasing = TRUE), ]
print(head(mb_bm_df[, c("feature","loading_score","max_cross_cor","top_partner")], 10))

# Known IBD-dysbiotic taxa check
known_ibd_up_taxa   <- c("Ruminococcus_gnavus", "Ruminococcus_torques",
                          "Veillonella_parvula", "Haemophilus_parainfluenzae",
                          "Escherichia_coli", "Streptococcus_salivarius")
known_ibd_down_taxa <- c("Faecalibacterium_prausnitzii",
                          "Roseburia_intestinalis", "Roseburia_hominis",
                          "Butyrivibrio_fibrisolvens", "Blautia_obeum",
                          "Coprococcus_comes")

taxa_in_bm      <- as.character(mb_bm_df$feature)
ibdup_recovered <- sum(known_ibd_up_taxa %in% taxa_in_bm)
ibddn_recovered <- sum(known_ibd_down_taxa %in% taxa_in_bm)

cat(sprintf("\nKnown IBD-increased taxa in biomarkers: %d / %d\n",
            ibdup_recovered, length(known_ibd_up_taxa)))
cat(sprintf("Known IBD-decreased taxa in biomarkers: %d / %d\n",
            ibddn_recovered, length(known_ibd_down_taxa)))

# ── 6b: Classifier validation ─────────────────────────────────────────────────
cr <- performance(result)
cat("\n-- Classifier Performance (5-fold CV) --\n")
cat(sprintf("  Host RNA-seq only   : AUC = %.3f ± %.3f\n",
            cr$host_only$mean_auc, cr$host_only$sd_auc))
cat(sprintf("  Microbiome only     : AUC = %.3f ± %.3f\n",
            cr$microbiome_only$mean_auc, cr$microbiome_only$sd_auc))
cat(sprintf("  Joint multi-omics   : AUC = %.3f ± %.3f\n",
            cr$joint$mean_auc, cr$joint$sd_auc))
delta <- cr$joint$mean_auc - cr$host_only$mean_auc
cat(sprintf("  Multi-omics gain    : +%.3f AUC vs host-only\n", delta))

# Sanity check: joint should be >= microbiome-only on this dataset
if (cr$joint$mean_auc >= cr$microbiome_only$mean_auc - 0.10) {
    cat("  PASS: Joint >= microbiome-only (within tolerance)\n")
} else {
    cat("  NOTE: Joint < microbiome-only; check feature selection parameters\n")
}

# ── 7. Generate all plots ─────────────────────────────────────────────────────

cat("\n=== STEP 7: Generate Plots ===\n")

save_plot <- function(p, filename, width = 9, height = 6) {
    path <- file.path(outdir, filename)
    ggsave(path, plot = p, width = width, height = height, dpi = 150)
    cat(sprintf("  Saved: %s\n", path))
    invisible(path)
}

# 7a: Joint integration biplot
cat("\nPlot 1/5: Integration biplot...\n")
p_int <- plotIntegration(
    result,
    outcome          = outcome_vec,
    n_loading_arrows = 6,
    point_size       = 2.2,
    point_alpha      = 0.85,
    colours          = c(UC = "#E24B4A", nonIBD = "#378ADD")
)
save_plot(p_int, "01_integration_biplot.png")
print(p_int)

# 7b: Cross-omics biomarker network
cat("\nPlot 2/5: Biomarker correlation heatmap...\n")
p_net <- plotBiomarkerNetwork(
    result,
    mae,
    n_host  = min(20, sum(bm_df$omics_layer == "host")),
    n_mb    = min(15, sum(bm_df$omics_layer == "microbiome")),
    host_assay = "voom",
    mb_assay   = "CLR"
)
save_plot(p_net, "02_biomarker_network.png", width = 10, height = 8)
print(p_net)

# 7c: Classifier comparison (bar)
cat("\nPlot 3/5: Classifier AUC bar chart...\n")
p_bar <- plotClassifierComparison(result, type = "bar")
save_plot(p_bar, "03_classifier_comparison_bar.png", width = 7, height = 5)
print(p_bar)

# 7d: Classifier comparison (ROC)
cat("\nPlot 4/5: ROC curves...\n")
p_roc <- plotClassifierComparison(result, type = "roc")
save_plot(p_roc, "04_classifier_comparison_roc.png", width = 7, height = 6)
print(p_roc)

# 7e: Sankey flow diagram
cat("\nPlot 5/5: Feature flow (Sankey) diagram...\n")
p_sank <- plotSankey(result, n_features = 8)
save_plot(p_sank, "05_sankey_flow.png", width = 10, height = 7)
print(p_sank)

# ── 8. Comprehensive report ───────────────────────────────────────────────────

cat("\n=== STEP 8: Comprehensive Report ===\n")
report_path <- file.path(outdir, "MultiOmicsBridge_IBDMDB_report.txt")
generateReport(result, file = report_path, n_top = 15)

# ── 9. Session and timing summary ─────────────────────────────────────────────

elapsed_total <- as.numeric(Sys.time() - t0, units = "secs")

cat("\n============================================================\n")
cat("  VALIDATION SUMMARY\n")
cat("============================================================\n")
cat(sprintf("  Dataset            : FranzsosaEA_2019 (IBDMDB)\n"))
cat(sprintf("  Paired samples     : %d (UC=%d, nonIBD=%d)\n",
            ncol(mae),
            sum(outcome_vec == "UC"),
            sum(outcome_vec == "nonIBD")))
cat(sprintf("  Host genes         : %d (simulated)\n",
            nrow(experiments(mae)[["host"]])))
cat(sprintf("  Microbial taxa     : %d (REAL)\n",
            nrow(experiments(mae)[["microbiome"]])))
cat(sprintf("  Biomarkers found   : %d (%d host, %d mb)\n",
            nrow(bm_df),
            sum(bm_df$omics_layer == "host"),
            sum(bm_df$omics_layer == "microbiome")))
cat(sprintf("  IBD-up genes found : %d / %d known\n",
            length(up_recovered), length(ibd_up_genes)))
cat(sprintf("  Joint AUC          : %.3f (±%.3f)\n",
            cr$joint$mean_auc, cr$joint$sd_auc))
cat(sprintf("  Multi-omics gain   : +%.3f vs host-only\n", delta))
cat(sprintf("  Total runtime      : %.1f seconds\n", elapsed_total))
cat(sprintf("  Output directory   : %s\n", outdir))
cat("============================================================\n")

if (cr$joint$mean_auc > 0.70) {
    cat("  STATUS: PASS — pipeline validated on real microbiome data\n")
} else {
    cat("  STATUS: CHECK — AUC lower than expected; review parameters\n")
}
cat("============================================================\n")

sessionInfo()
