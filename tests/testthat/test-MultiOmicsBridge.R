library(testthat)
library(MultiOmicsBridge)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(S4Vectors)

# ── Shared test helpers ────────────────────────────────────────────────────────

#' Minimal host count matrix (200 genes x 20 samples)
make_host_counts <- function(seed = 42, n_genes = 200, n_samples = 20) {
  set.seed(seed)
  counts <- matrix(rpois(n_genes * n_samples, lambda = 150L),
    nrow = n_genes, ncol = n_samples
  )
  rownames(counts) <- paste0("Gene", seq_len(n_genes))
  colnames(counts) <- paste0("Sample", seq_len(n_samples))
  # Inject modest DE signal for first 20 genes in treatment group
  counts[seq_len(20), seq(11, n_samples)] <-
    counts[seq_len(20), seq(11, n_samples)] * 4L
  counts
}

#' Minimal microbiome count matrix (60 taxa x 20 samples)
make_mb_counts <- function(seed = 42, n_taxa = 60, n_samples = 20) {
  set.seed(seed)
  counts <- matrix(rpois(n_taxa * n_samples, lambda = 40L),
    nrow = n_taxa, ncol = n_samples
  )
  rownames(counts) <- paste0("Taxon", seq_len(n_taxa))
  colnames(counts) <- paste0("Sample", seq_len(n_samples))
  # Inject signal in first 5 taxa for treatment group
  counts[seq_len(5), seq(11, n_samples)] <-
    counts[seq_len(5), seq(11, n_samples)] * 3L
  counts
}

make_col_data <- function(n_samples = 20) {
  data.frame(
    condition = rep(c("ctrl", "treat"), each = n_samples / 2),
    age       = sample(30:60, n_samples, replace = TRUE),
    row.names = paste0("Sample", seq_len(n_samples))
  )
}

# Full pipeline helper (used across many test blocks)
make_full_result <- function(seed = 42, cv_folds = 3L) {
  host_counts <- make_host_counts(seed)
  mb_counts <- make_mb_counts(seed)
  col_data <- make_col_data()

  host_se <- loadHostData(host_counts, col_data = col_data)
  mb_se <- loadMicrobiomeData(mb_counts)
  mae <- matchSamples(host_se, mb_se)
  outcome <- rep(c("ctrl", "treat"), each = 10)

  MultiOmicsBridgeAnalysis(
    mae,
    outcome,
    n_components    = 2L,
    n_features_host = 15L,
    n_features_mb   = 10L,
    n_biomarkers    = 10L,
    cv_folds        = cv_folds,
    seed            = seed
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# MOBResult S4 class and accessors
# ══════════════════════════════════════════════════════════════════════════════

test_that("MOBResult constructor works", {
  scores <- matrix(rnorm(20),
    nrow = 10, ncol = 2,
    dimnames = list(
      paste0("S", 1:10),
      c("Comp1", "Comp2")
    )
  )
  fl <- list(
    host = matrix(rnorm(10),
      nrow = 5, ncol = 2,
      dimnames = list(paste0("G", 1:5), c("Comp1", "Comp2"))
    ),
    microbiome = matrix(rnorm(6),
      nrow = 3, ncol = 2,
      dimnames = list(paste0("T", 1:3), c("Comp1", "Comp2"))
    )
  )
  bm <- DataFrame(
    feature = c("G1", "T1"),
    omics_layer = c("host", "microbiome"),
    loading_score = c(0.8, 0.6),
    rank = c(1L, 2L),
    component = c(1L, 1L)
  )
  obj <- MOBResult(scores, fl, bm, list(),
    params = list(
      integration_method = "DIABLO",
      outcome_levels = c("ctrl", "treat")
    )
  )
  expect_s4_class(obj, "MOBResult")
})

test_that("MOBResult accessors return correct slot types", {
  scores <- matrix(rnorm(20), 10, 2)
  fl <- list(
    host = matrix(rnorm(10), 5, 2),
    microbiome = matrix(rnorm(6), 3, 2)
  )
  bm <- DataFrame(
    feature = c("G1", "T1"),
    omics_layer = c("host", "microbiome"),
    loading_score = c(0.8, 0.6),
    rank = 1:2, component = c(1L, 1L)
  )
  cr <- list(
    host_only = list(mean_auc = 0.8),
    microbiome_only = list(mean_auc = 0.7),
    joint = list(mean_auc = 0.9)
  )
  obj <- MOBResult(scores, fl, bm, cr)

  expect_true(is.matrix(integrationScores(obj)))
  expect_type(featureLoadings(obj), "list")
  expect_s4_class(biomarkers(obj), "DataFrame")
  expect_type(performance(obj), "list")
})

test_that("show method for MOBResult prints without error", {
  scores <- matrix(rnorm(20), 10, 2)
  obj <- MOBResult(scores, list(), DataFrame(), list(),
    params = list(
      outcome_levels = c("ctrl", "treat"),
      integration_method = "DIABLO"
    )
  )
  expect_output(show(obj), "MOBResult")
})

# ══════════════════════════════════════════════════════════════════════════════
# loadHostData
# ══════════════════════════════════════════════════════════════════════════════

test_that("loadHostData returns a SummarizedExperiment", {
  counts <- make_host_counts()
  host_se <- loadHostData(counts)
  expect_s4_class(host_se, "SummarizedExperiment")
})

test_that("loadHostData has 'counts' and 'voom' assays", {
  counts <- make_host_counts()
  host_se <- loadHostData(counts)
  expect_true("counts" %in% assayNames(host_se))
  expect_true("voom" %in% assayNames(host_se))
})

test_that("loadHostData preserves gene names", {
  counts <- make_host_counts()
  host_se <- loadHostData(counts)
  expect_equal(rownames(host_se), rownames(counts))
})

test_that("loadHostData includes colData when provided", {
  counts <- make_host_counts()
  col_data <- make_col_data()
  host_se <- loadHostData(counts, col_data = col_data)
  expect_true("condition" %in% names(colData(host_se)))
})

test_that("loadHostData errors on negative counts", {
  counts <- make_host_counts()
  counts[1, 1] <- -1L
  expect_error(loadHostData(counts), regexp = "negative")
})

test_that("loadHostData errors on missing rownames", {
  counts <- make_host_counts()
  rownames(counts) <- NULL
  expect_error(loadHostData(counts), regexp = "rownames")
})

test_that("loadHostData accepts SummarizedExperiment input", {
  counts <- make_host_counts()
  se_in <- SummarizedExperiment(assays = list(counts = counts))
  host_se <- loadHostData(se_in)
  expect_s4_class(host_se, "SummarizedExperiment")
})

test_that("loadHostData min_count filter removes all-zero genes", {
  counts <- make_host_counts()
  counts[1, ] <- 0L # set first gene to all-zero
  host_se <- loadHostData(counts, min_count = 1)
  expect_false("Gene1" %in% rownames(host_se))
})

# ══════════════════════════════════════════════════════════════════════════════
# loadMicrobiomeData
# ══════════════════════════════════════════════════════════════════════════════

test_that("loadMicrobiomeData returns a SummarizedExperiment", {
  counts <- make_mb_counts()
  mb_se <- loadMicrobiomeData(counts)
  expect_s4_class(mb_se, "SummarizedExperiment")
})

test_that("loadMicrobiomeData with CLR has 'counts' and 'CLR' assays", {
  counts <- make_mb_counts()
  mb_se <- loadMicrobiomeData(counts, normalization = "CLR")
  expect_true("counts" %in% assayNames(mb_se))
  expect_true("CLR" %in% assayNames(mb_se))
})

test_that("loadMicrobiomeData with TSS has 'TSS' assay", {
  counts <- make_mb_counts()
  mb_se <- loadMicrobiomeData(counts, normalization = "TSS")
  expect_true("TSS" %in% assayNames(mb_se))
})

test_that("CLR transformation produces zero-mean columns", {
  counts <- make_mb_counts()
  mb_se <- loadMicrobiomeData(counts,
    normalization = "CLR",
    min_prevalence = 0
  )
  clr_mat <- assay(mb_se, "CLR")
  # Each column (sample) should have mean ≈ 0
  col_means <- abs(colMeans(clr_mat))
  expect_true(all(col_means < 1e-9))
})

test_that("loadMicrobiomeData removes low-prevalence taxa", {
  counts <- make_mb_counts()
  counts[1, ] <- 0L # zero out first taxon
  mb_se <- loadMicrobiomeData(counts, min_prevalence = 0.05)
  expect_false("Taxon1" %in% rownames(mb_se))
})

test_that("loadMicrobiomeData errors on negative values", {
  counts <- make_mb_counts()
  counts[1, 1] <- -1L
  expect_error(loadMicrobiomeData(counts), regexp = "negative")
})

test_that("loadMicrobiomeData errors on bad min_prevalence", {
  counts <- make_mb_counts()
  expect_error(loadMicrobiomeData(counts, min_prevalence = 1.5),
    regexp = "\\[0, 1\\]"
  )
})

# ══════════════════════════════════════════════════════════════════════════════
# matchSamples
# ══════════════════════════════════════════════════════════════════════════════

test_that("matchSamples returns a MultiAssayExperiment", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  expect_s4_class(mae, "MultiAssayExperiment")
})

test_that("matchSamples retains correct common samples", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  # All 20 samples should match
  expect_equal(nrow(colData(mae)), 20L)
})

test_that("matchSamples handles partially overlapping samples", {
  # Host: Sample1-20, Microbiome: Sample1-15 (5 missing)
  host_counts <- make_host_counts()
  mb_counts <- make_mb_counts()[, 1:15]

  host_se <- loadHostData(host_counts)
  mb_se <- loadMicrobiomeData(mb_counts)
  mae <- matchSamples(host_se, mb_se, min_paired = 5)
  expect_equal(nrow(colData(mae)), 15L)
})

test_that("matchSamples errors when no common samples", {
  host_counts <- make_host_counts()
  mb_counts <- make_mb_counts()
  colnames(mb_counts) <- paste0("DifferentSample", 1:20)

  host_se <- loadHostData(host_counts)
  mb_se <- loadMicrobiomeData(mb_counts)
  expect_error(matchSamples(host_se, mb_se), regexp = "No common")
})

test_that("matchSamples errors below min_paired threshold", {
  host_counts <- make_host_counts()
  mb_counts <- make_mb_counts()[, 1:3] # only 3 shared samples

  host_se <- loadHostData(host_counts)
  mb_se <- loadMicrobiomeData(mb_counts)
  expect_error(
    matchSamples(host_se, mb_se, min_paired = 10),
    regexp = "Only 3 paired"
  )
})

test_that("matchSamples MAE has 'host' and 'microbiome' experiments", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  expect_true("host" %in% names(experiments(mae)))
  expect_true("microbiome" %in% names(experiments(mae)))
})

# ══════════════════════════════════════════════════════════════════════════════
# jointDimReduction
# ══════════════════════════════════════════════════════════════════════════════

test_that("jointDimReduction returns a named list", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  outcome <- rep(c("ctrl", "treat"), each = 10)
  dr <- jointDimReduction(mae, outcome,
    n_components = 2,
    n_features_host = 10,
    n_features_mb = 8
  )
  expect_type(dr, "list")
  expect_true(all(c(
    "scores", "host_loadings", "mb_loadings",
    "explained_variance"
  ) %in% names(dr)))
})

test_that("jointDimReduction scores have correct dimensions", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  outcome <- rep(c("ctrl", "treat"), each = 10)
  dr <- jointDimReduction(mae, outcome,
    n_components = 2,
    n_features_host = 10,
    n_features_mb = 8
  )
  expect_equal(nrow(dr$scores), 20L)
  expect_equal(ncol(dr$scores), 2L)
})

test_that("jointDimReduction errors on non-MAE input", {
  expect_error(
    jointDimReduction(list(), outcome = "ctrl"),
    regexp = "MultiAssayExperiment"
  )
})

test_that("jointDimReduction errors when outcome length mismatches samples", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  expect_error(
    jointDimReduction(mae,
      outcome = c("ctrl", "treat"),
      n_features_host = 10, n_features_mb = 8
    ),
    regexp = "Length of"
  )
})

test_that("jointDimReduction errors with single outcome level", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  expect_error(
    jointDimReduction(mae,
      outcome = rep("ctrl", 20),
      n_features_host = 10, n_features_mb = 8
    ),
    regexp = "at least 2 levels"
  )
})

# ══════════════════════════════════════════════════════════════════════════════
# biomarkerDiscovery
# ══════════════════════════════════════════════════════════════════════════════

test_that("biomarkerDiscovery returns a DataFrame", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  outcome <- rep(c("ctrl", "treat"), each = 10)
  dr <- jointDimReduction(mae, outcome,
    n_components = 2,
    n_features_host = 10, n_features_mb = 8
  )
  bm <- biomarkerDiscovery(mae, dr, n_biomarkers = 8)
  expect_s4_class(bm, "DataFrame")
})

test_that("biomarkerDiscovery output has required columns", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  outcome <- rep(c("ctrl", "treat"), each = 10)
  dr <- jointDimReduction(mae, outcome,
    n_components = 2,
    n_features_host = 10, n_features_mb = 8
  )
  bm <- biomarkerDiscovery(mae, dr, n_biomarkers = 8)
  expect_true(all(c(
    "feature", "omics_layer", "loading_score",
    "rank", "component"
  ) %in% names(bm)))
})

test_that("biomarkerDiscovery includes both host and microbiome features", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  outcome <- rep(c("ctrl", "treat"), each = 10)
  dr <- jointDimReduction(mae, outcome,
    n_components = 2,
    n_features_host = 10, n_features_mb = 8
  )
  bm <- biomarkerDiscovery(mae, dr, n_biomarkers = 8)
  layers <- unique(as.character(bm[["omics_layer"]]))
  expect_true("host" %in% layers)
  expect_true("microbiome" %in% layers)
})

test_that("biomarkerDiscovery errors on invalid dr_result", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  expect_error(biomarkerDiscovery(mae, list(bad = 1)),
    regexp = "output of jointDimReduction"
  )
})

# ══════════════════════════════════════════════════════════════════════════════
# diagnosticClassifier
# ══════════════════════════════════════════════════════════════════════════════

test_that("diagnosticClassifier returns a named list", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  outcome <- rep(c("ctrl", "treat"), each = 10)
  dr <- jointDimReduction(mae, outcome,
    n_components = 2,
    n_features_host = 10, n_features_mb = 8
  )
  bm <- biomarkerDiscovery(mae, dr, n_biomarkers = 8)
  clf <- diagnosticClassifier(mae, outcome,
    biomarker_table = bm,
    cv_folds = 3L
  )
  expect_type(clf, "list")
  expect_true(all(c("host_only", "microbiome_only", "joint") %in% names(clf)))
})

test_that("diagnosticClassifier AUC values are in [0, 1]", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  outcome <- rep(c("ctrl", "treat"), each = 10)
  dr <- jointDimReduction(mae, outcome,
    n_components = 2,
    n_features_host = 10, n_features_mb = 8
  )
  bm <- biomarkerDiscovery(mae, dr, n_biomarkers = 8)
  clf <- diagnosticClassifier(mae, outcome,
    biomarker_table = bm,
    cv_folds = 3L
  )
  for (nm in c("host_only", "microbiome_only", "joint")) {
    auc <- clf[[nm]]$mean_auc
    expect_gte(auc, 0)
    expect_lte(auc, 1)
  }
})

test_that("joint AUC is at least as high as microbiome-only on signal data", {
  # On data with injected multi-omics signal, joint should >= mb-only
  clf <- make_full_result()
  cr <- performance(clf)
  expect_gte(cr$joint$mean_auc, cr$microbiome_only$mean_auc - 0.05)
})

test_that("diagnosticClassifier errors on non-binary outcome", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  outcome <- rep(c("a", "b", "c"), length.out = 20)
  expect_error(
    diagnosticClassifier(mae, outcome, cv_folds = 3L),
    regexp = "exactly 2 levels"
  )
})

test_that("diagnosticClassifier n_features slot is filled correctly", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  outcome <- rep(c("ctrl", "treat"), each = 10)
  dr <- jointDimReduction(mae, outcome,
    n_components = 2,
    n_features_host = 10, n_features_mb = 8
  )
  bm <- biomarkerDiscovery(mae, dr, n_biomarkers = 8)
  clf <- diagnosticClassifier(mae, outcome,
    biomarker_table = bm,
    cv_folds = 3L
  )
  expect_true("host" %in% names(clf$n_features))
  expect_true(all(clf$n_features > 0))
})

# ══════════════════════════════════════════════════════════════════════════════
# MultiOmicsBridgeAnalysis (main wrapper)
# ══════════════════════════════════════════════════════════════════════════════

test_that("MultiOmicsBridgeAnalysis returns a MOBResult", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  outcome <- rep(c("ctrl", "treat"), each = 10)
  result <- MultiOmicsBridgeAnalysis(mae, outcome,
    n_components = 2,
    n_features_host = 10,
    n_features_mb = 8,
    n_biomarkers = 8,
    cv_folds = 3L
  )
  expect_s4_class(result, "MOBResult")
})

test_that("MultiOmicsBridgeAnalysis result has non-empty biomarker table", {
  result <- make_full_result()
  expect_gt(nrow(biomarkers(result)), 0L)
})

test_that("MultiOmicsBridgeAnalysis result has three classifier entries", {
  result <- make_full_result()
  cr <- performance(result)
  expect_true(all(c("host_only", "microbiome_only", "joint") %in% names(cr)))
})

test_that("MultiOmicsBridgeAnalysis params are stored correctly", {
  result <- make_full_result()
  p <- result@params
  expect_equal(p$integration_method, "DIABLO")
  expect_equal(p$cv_folds, 3L)
  expect_equal(p$outcome_levels, c("ctrl", "treat"))
})

test_that("MultiOmicsBridgeAnalysis errors on non-MAE input", {
  expect_error(
    MultiOmicsBridgeAnalysis(list(), outcome = rep("ctrl", 10)),
    regexp = "MultiAssayExperiment"
  )
})

test_that("MultiOmicsBridgeAnalysis errors on mismatched outcome length", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  expect_error(
    MultiOmicsBridgeAnalysis(mae,
      outcome = c("ctrl", "treat"),
      n_features_host = 5, n_features_mb = 5,
      cv_folds = 3L
    ),
    regexp = "must match"
  )
})

test_that("MultiOmicsBridgeAnalysis errors on single-level outcome", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  expect_error(
    MultiOmicsBridgeAnalysis(mae,
      outcome = rep("ctrl", 20),
      n_features_host = 5, n_features_mb = 5,
      cv_folds = 3L
    ),
    regexp = "exactly 2 levels"
  )
})

# ══════════════════════════════════════════════════════════════════════════════
# Visualization functions
# ══════════════════════════════════════════════════════════════════════════════

test_that("plotIntegration returns a ggplot object", {
  result <- make_full_result()
  outcome <- rep(c("ctrl", "treat"), each = 10)
  p <- plotIntegration(result, outcome = outcome)
  expect_s3_class(p, "ggplot")
})

test_that("plotIntegration errors on non-MOBResult input", {
  expect_error(plotIntegration(list()), regexp = "MOBResult")
})

test_that("plotBiomarkerNetwork returns a ggplot object", {
  host_se <- loadHostData(make_host_counts())
  mb_se <- loadMicrobiomeData(make_mb_counts())
  mae <- matchSamples(host_se, mb_se)
  result <- make_full_result()
  p <- plotBiomarkerNetwork(result, mae, n_host = 5, n_mb = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plotBiomarkerNetwork errors on non-MOBResult input", {
  mae <- matchSamples(
    loadHostData(make_host_counts()),
    loadMicrobiomeData(make_mb_counts())
  )
  expect_error(plotBiomarkerNetwork(list(), mae), regexp = "MOBResult")
})

test_that("plotClassifierComparison returns a ggplot (bar type)", {
  result <- make_full_result()
  p <- plotClassifierComparison(result, type = "bar")
  expect_s3_class(p, "ggplot")
})

test_that("plotClassifierComparison returns a ggplot (roc type)", {
  result <- make_full_result()
  p <- plotClassifierComparison(result, type = "roc")
  expect_s3_class(p, "ggplot")
})

test_that("plotClassifierComparison errors on non-MOBResult input", {
  expect_error(plotClassifierComparison(list()), regexp = "MOBResult")
})

test_that("plotSankey returns a ggplot object", {
  result <- make_full_result()
  p <- plotSankey(result, n_features = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plotSankey errors on non-MOBResult input", {
  expect_error(plotSankey(list()), regexp = "MOBResult")
})

# ══════════════════════════════════════════════════════════════════════════════
# generateReport
# ══════════════════════════════════════════════════════════════════════════════

test_that("generateReport prints without error", {
  result <- make_full_result()
  expect_message(generateReport(result, n_top = 5), "MultiOmicsBridge")
})

test_that("generateReport saves file when path is given", {
  result <- make_full_result()
  tmp_file <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp_file))
  capture.output(suppressMessages(generateReport(result, file = tmp_file, n_top = 3)))
  expect_true(file.exists(tmp_file))
  content <- readLines(tmp_file)
  expect_true(any(grepl("MultiOmicsBridge", content)))
})

test_that("generateReport returns a named list invisibly", {
  result <- make_full_result()
  out <- NULL
  capture.output(out <- generateReport(result, n_top = 3))
  expect_type(out, "list")
  expect_true("header" %in% names(out))
})

test_that("generateReport errors on non-MOBResult input", {
  expect_error(generateReport(list()), regexp = "MOBResult")
})
