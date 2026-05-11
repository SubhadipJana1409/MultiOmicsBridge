#' @title Full Multi-Omics Integration Analysis Pipeline
#'
#' @description
#' A one-call wrapper that executes the complete MultiOmicsBridge analysis
#' pipeline: joint dimensionality reduction, multi-omics biomarker
#' discovery, and integrated diagnostic classification. The function
#' accepts a \code{MultiAssayExperiment} (from \code{\link{matchSamples}})
#' and returns a \code{\link{MOBResult}} S4 object containing all results.
#'
#' @details
#' Internally, \code{MultiOmicsBridgeAnalysis} calls:
#' \enumerate{
#'   \item \code{\link{jointDimReduction}} — DIABLO sparse multi-block
#'     PLS-DA.
#'   \item \code{\link{biomarkerDiscovery}} — loading-based ranking and
#'     cross-omics correlation annotation.
#'   \item \code{\link{diagnosticClassifier}} — host-only, microbiome-only,
#'     and joint Random Forest classifiers with cross-validation.
#' }
#' Each step can also be called independently for fine-grained control.
#'
#' @param mae A \code{MultiAssayExperiment} with \code{"host"} and
#'   \code{"microbiome"} experiments, produced by \code{\link{matchSamples}}.
#' @param outcome A \code{character} or \code{factor} vector of outcome
#'   labels, one per sample (ordered as \code{colnames(mae)}).
#'   Must have exactly 2 levels for classification.
#' @param n_components An \code{integer(1)} number of DIABLO latent
#'   components to extract. Default: \code{2}.
#' @param n_features_host An \code{integer(1)} number of host genes
#'   to select per DIABLO component. Default: \code{50}.
#' @param n_features_mb An \code{integer(1)} number of microbial taxa
#'   to select per DIABLO component. Default: \code{20}.
#' @param n_biomarkers An \code{integer(1)} total number of biomarkers
#'   to report per omics layer in the ranked table. Default: \code{50}.
#' @param cv_folds An \code{integer(1)} number of cross-validation folds
#'   for classifier evaluation. Default: \code{5}.
#' @param host_assay A \code{character(1)} name of the normalized host
#'   assay in the MAE. Default: \code{"voom"}.
#' @param mb_assay A \code{character(1)} name of the normalized microbiome
#'   assay in the MAE. Default: \code{"CLR"}.
#' @param design_off_diag A \code{numeric(1)} off-diagonal value for the
#'   DIABLO design matrix (see \code{\link{jointDimReduction}}).
#'   Default: \code{0.1}.
#' @param seed An \code{integer(1)} random seed. Default: \code{42}.
#' @param BPPARAM A \code{BiocParallelParam} object. Default:
#'   \code{SerialParam()}.
#'
#' @return A \code{\link{MOBResult}} object with:
#'   \describe{
#'     \item{\code{integrationScores(result)}}{Matrix of DIABLO sample
#'       scores (samples x components).}
#'     \item{\code{featureLoadings(result)}}{Named list of host and
#'       microbiome loading matrices.}
#'     \item{\code{biomarkers(result)}}{Ranked multi-omics biomarker
#'       \code{DataFrame}.}
#'     \item{\code{performance(result)}}{Classifier performance list with
#'       AUC-ROC for host-only, microbiome-only, and joint models.}
#'   }
#'
#' @examples
#' set.seed(42)
#' host_counts <- matrix(rpois(500 * 20, 100), nrow = 500, ncol = 20,
#'     dimnames = list(paste0("Gene", 1:500), paste0("S", 1:20)))
#' host_counts[1:20, 11:20] <- host_counts[1:20, 11:20] * 5L
#' mb_counts <- matrix(rpois(60 * 20, 40), nrow = 60, ncol = 20,
#'     dimnames = list(paste0("Taxon", 1:60), paste0("S", 1:20)))
#'
#' host_se <- loadHostData(host_counts)
#' mb_se   <- loadMicrobiomeData(mb_counts)
#' mae     <- matchSamples(host_se, mb_se)
#' outcome <- rep(c("ctrl", "treat"), each = 10)
#'
#' result <- MultiOmicsBridgeAnalysis(mae, outcome,
#'                                     n_components    = 2,
#'                                     n_features_host = 20,
#'                                     n_features_mb   = 10,
#'                                     n_biomarkers    = 15,
#'                                     cv_folds        = 3)
#' result
#' head(as.data.frame(biomarkers(result)))
#'
#' @seealso \code{\link{matchSamples}}, \code{\link{jointDimReduction}},
#'   \code{\link{biomarkerDiscovery}}, \code{\link{diagnosticClassifier}},
#'   \code{\link{plotIntegration}}, \code{\link{plotClassifierComparison}}
#'
#' @importFrom MultiAssayExperiment experiments colData
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is
#' @export
MultiOmicsBridgeAnalysis <- function(mae,
                                      outcome,
                                      n_components     = 2L,
                                      n_features_host  = 50L,
                                      n_features_mb    = 20L,
                                      n_biomarkers     = 50L,
                                      cv_folds         = 5L,
                                      host_assay       = "voom",
                                      mb_assay         = "CLR",
                                      design_off_diag  = 0.1,
                                      seed             = 42L,
                                      BPPARAM          = SerialParam()) {
    # ── Validation ────────────────────────────────────────────────────────────
    if (!is(mae, "MultiAssayExperiment"))
        stop("'mae' must be a MultiAssayExperiment. ",
             "Run matchSamples() first.")

    exp_names <- names(MultiAssayExperiment::experiments(mae))
    if (!all(c("host", "microbiome") %in% exp_names))
        stop("'mae' must contain 'host' and 'microbiome' experiments.")

    n_samples <- nrow(MultiAssayExperiment::colData(mae))
    if (length(outcome) != n_samples)
        stop("'outcome' length (", length(outcome),
             ") must match number of samples (", n_samples, ").")

    outcome_f <- factor(outcome)
    cond_levels <- levels(outcome_f)
    if (nlevels(outcome_f) != 2L)
        stop("'outcome' must have exactly 2 levels. Found: ",
             paste(cond_levels, collapse = ", "))

    message("=== MultiOmicsBridgeAnalysis ===")
    message(sprintf("Samples: %d | Outcome: %s vs %s",
                    n_samples, cond_levels[1L], cond_levels[2L]))

    # ── Module 2: Joint Dimensionality Reduction ──────────────────────────────
    message("\n-- Module 2: Joint Dimensionality Reduction (DIABLO) --")
    dr_result <- jointDimReduction(
        mae             = mae,
        outcome         = outcome_f,
        n_components    = n_components,
        n_features_host = n_features_host,
        n_features_mb   = n_features_mb,
        design_off_diag = design_off_diag,
        host_assay      = host_assay,
        mb_assay        = mb_assay
    )

    # ── Module 3: Biomarker Discovery ─────────────────────────────────────────
    message("\n-- Module 3: Biomarker Discovery --")
    bm_table <- biomarkerDiscovery(
        mae          = mae,
        dr_result    = dr_result,
        n_biomarkers = n_biomarkers,
        host_assay   = host_assay,
        mb_assay     = mb_assay
    )

    # ── Module 4: Diagnostic Classification ──────────────────────────────────
    message("\n-- Module 4: Diagnostic Classification --")
    clf_results <- diagnosticClassifier(
        mae             = mae,
        outcome         = outcome_f,
        biomarker_table = bm_table,
        cv_folds        = cv_folds,
        seed            = seed,
        host_assay      = host_assay,
        mb_assay        = mb_assay
    )

    # ── Assemble MOBResult ────────────────────────────────────────────────────
    message("\n=== Analysis complete. Assembling MOBResult ===")

    # Pad scores to n_samples x n_components matrix
    scores <- if (!is.null(dr_result$scores)) {
        dr_result$scores
    } else {
        matrix(nrow = n_samples, ncol = n_components,
               dimnames = list(colnames(mae),
                               paste0("Comp", seq_len(n_components))))
    }

    MOBResult(
        integratedScores  = scores,
        featureLoadings   = list(
            host       = dr_result$host_loadings,
            microbiome = dr_result$mb_loadings
        ),
        biomarkerTable    = bm_table,
        classifierResults = clf_results,
        params = list(
            n_components     = n_components,
            n_features_host  = n_features_host,
            n_features_mb    = n_features_mb,
            n_biomarkers     = n_biomarkers,
            cv_folds         = cv_folds,
            host_assay       = host_assay,
            mb_assay         = mb_assay,
            integration_method = "DIABLO",
            design_off_diag  = design_off_diag,
            outcome_levels   = cond_levels,
            seed             = seed
        )
    )
}
