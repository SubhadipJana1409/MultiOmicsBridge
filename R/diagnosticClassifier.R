#' @title Integrated Diagnostic Classifier with Multi-Omics Comparison
#'
#' @description
#' Trains and evaluates three Random Forest diagnostic classifiers using
#' cross-validation: a host-only model, a microbiome-only model, and a
#' joint multi-omics model. By comparing AUC-ROC across all three
#' configurations, \code{diagnosticClassifier} quantifies the added
#' diagnostic value of combining both data types.
#'
#' @details
#' Each classifier is trained using \code{ranger::ranger} (fast C++
#' Random Forest) with 500 trees and stratified k-fold cross-validation.
#' Features are drawn from the top biomarkers identified by
#' \code{\link{biomarkerDiscovery}} (or all available features if
#' \code{biomarker_table} is \code{NULL}).
#'
#' The cross-validation procedure:
#' \enumerate{
#'   \item Split samples into \code{cv_folds} stratified folds (each
#'     fold preserves the outcome class ratio).
#'   \item For each fold, train on the remaining folds, predict on the
#'     held-out fold.
#'   \item Compute AUC-ROC on held-out predictions.
#'   \item Report mean ± SD AUC across folds.
#' }
#'
#' @param mae A \code{MultiAssayExperiment} from \code{\link{matchSamples}}.
#' @param outcome A \code{character} or \code{factor} vector of outcome
#'   labels, one per sample. Must have exactly 2 levels.
#' @param biomarker_table An optional \code{DataFrame} from
#'   \code{\link{biomarkerDiscovery}}. If provided, only the listed
#'   features are used as classifier inputs. If \code{NULL}, all features
#'   in the matched MAE are used. Default: \code{NULL}.
#' @param cv_folds An \code{integer(1)} number of cross-validation folds.
#'   Default: \code{5}.
#' @param n_trees An \code{integer(1)} number of trees in each Random
#'   Forest. Default: \code{500}.
#' @param seed An \code{integer(1)} random seed for reproducibility.
#'   Default: \code{42}.
#' @param host_assay A \code{character(1)} assay in the host SE. Default:
#'   \code{"voom"}.
#' @param mb_assay A \code{character(1)} assay in the microbiome SE.
#'   Default: \code{"CLR"}.
#'
#' @return A named \code{list} with elements:
#'   \describe{
#'     \item{\code{host_only}}{List with \code{mean_auc}, \code{sd_auc},
#'       \code{fold_auc}, \code{roc_data} for the host-only classifier.}
#'     \item{\code{microbiome_only}}{Same structure for the
#'       microbiome-only classifier.}
#'     \item{\code{joint}}{Same structure for the joint classifier.}
#'     \item{\code{n_features}}{Named integer: \code{host},
#'       \code{microbiome}, \code{joint} feature counts.}
#'     \item{\code{cv_folds}}{Number of CV folds used.}
#'     \item{\code{outcome_levels}}{Outcome levels.}
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
#' dr_res  <- jointDimReduction(mae, outcome, n_components = 2,
#'                               n_features_host = 30, n_features_mb = 15)
#' bm      <- biomarkerDiscovery(mae, dr_res, n_biomarkers = 20)
#'
#' clf_res <- diagnosticClassifier(mae, outcome,
#'                                  biomarker_table = bm, cv_folds = 3)
#' clf_res$host_only$mean_auc
#' clf_res$joint$mean_auc
#'
#' @seealso \code{\link{biomarkerDiscovery}}, \code{\link{plotClassifierComparison}},
#'   \code{\link{MultiOmicsBridgeAnalysis}}
#'
#' @importFrom MultiAssayExperiment experiments colData
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom methods is
#' @importFrom stats var
#' @export
diagnosticClassifier <- function(mae,
                                  outcome,
                                  biomarker_table = NULL,
                                  cv_folds        = 5L,
                                  n_trees         = 500L,
                                  seed            = 42L,
                                  host_assay      = "voom",
                                  mb_assay        = "CLR") {
    # ── Validation ────────────────────────────────────────────────────────────
    if (!is(mae, "MultiAssayExperiment"))
        stop("'mae' must be a MultiAssayExperiment.")

    n_samples <- nrow(MultiAssayExperiment::colData(mae))
    if (length(outcome) != n_samples)
        stop("Length of 'outcome' (", length(outcome), ") must equal ",
             "number of MAE samples (", n_samples, ").")

    outcome_f <- factor(outcome)
    if (nlevels(outcome_f) != 2L)
        stop("'outcome' must have exactly 2 levels. Found: ",
             paste(levels(outcome_f), collapse = ", "))

    if (!is.numeric(cv_folds) || cv_folds < 2L)
        stop("'cv_folds' must be at least 2.")
    if (n_samples < cv_folds * 2L)
        stop("Too few samples (", n_samples, ") for ", cv_folds,
             "-fold CV. Reduce 'cv_folds'.")

    # ── Extract assay matrices ────────────────────────────────────────────────
    host_se <- MultiAssayExperiment::experiments(mae)[["host"]]
    mb_se   <- MultiAssayExperiment::experiments(mae)[["microbiome"]]

    host_assay_use <- if (host_assay %in% SummarizedExperiment::assayNames(host_se))
        host_assay else SummarizedExperiment::assayNames(host_se)[1L]
    mb_assay_use   <- if (mb_assay %in% SummarizedExperiment::assayNames(mb_se))
        mb_assay else SummarizedExperiment::assayNames(mb_se)[1L]

    host_mat <- t(SummarizedExperiment::assay(host_se, host_assay_use))
    mb_mat   <- t(SummarizedExperiment::assay(mb_se,   mb_assay_use))

    # ── Filter to biomarker features if provided ─────────────────────────────
    if (!is.null(biomarker_table) && nrow(biomarker_table) > 0L) {
        host_bm <- as.character(
            biomarker_table[["feature"]][
                biomarker_table[["omics_layer"]] == "host"])
        mb_bm   <- as.character(
            biomarker_table[["feature"]][
                biomarker_table[["omics_layer"]] == "microbiome"])

        host_bm_use <- intersect(host_bm, colnames(host_mat))
        mb_bm_use   <- intersect(mb_bm,   colnames(mb_mat))

        if (length(host_bm_use) > 0L)
            host_mat <- host_mat[, host_bm_use, drop = FALSE]
        if (length(mb_bm_use) > 0L)
            mb_mat   <- mb_mat[, mb_bm_use,     drop = FALSE]
    }

    # Remove zero-variance features
    host_var <- apply(host_mat, 2L, var, na.rm = TRUE)
    mb_var   <- apply(mb_mat,   2L, var, na.rm = TRUE)
    host_mat <- host_mat[, host_var > 0, drop = FALSE]
    mb_mat   <- mb_mat[, mb_var > 0,     drop = FALSE]

    n_host <- ncol(host_mat)
    n_mb   <- ncol(mb_mat)

    if (n_host == 0L)
        stop("No host features with non-zero variance remain.")
    if (n_mb == 0L)
        stop("No microbiome features with non-zero variance remain.")

    message(sprintf(
        "diagnosticClassifier: %d host | %d microbiome | %d-fold CV",
        n_host, n_mb, cv_folds
    ))

    # Joint feature matrix (column-bind)
    joint_mat <- cbind(
        host_mat,
        mb_mat
    )
    colnames(joint_mat) <- c(
        paste0("host___", colnames(host_mat)),
        paste0("mb___",   colnames(mb_mat))
    )

    # ── Stratified fold assignments ───────────────────────────────────────────
    fold_ids <- .stratifiedFolds(outcome_f, k = cv_folds, seed = seed)

    # ── Train all three classifiers ───────────────────────────────────────────
    message("  Training host-only classifier...")
    host_res <- .cvRandomForest(host_mat, outcome_f, fold_ids, seed = seed)

    message("  Training microbiome-only classifier...")
    mb_res <- .cvRandomForest(mb_mat, outcome_f, fold_ids, seed = seed)

    message("  Training joint classifier...")
    joint_res <- .cvRandomForest(joint_mat, outcome_f, fold_ids, seed = seed)

    message(sprintf(
        "AUC: host=%.3f | microbiome=%.3f | joint=%.3f",
        host_res$mean_auc, mb_res$mean_auc, joint_res$mean_auc
    ))

    list(
        host_only       = host_res,
        microbiome_only = mb_res,
        joint           = joint_res,
        n_features      = c(host       = n_host,
                            microbiome = n_mb,
                            joint      = n_host + n_mb),
        cv_folds        = cv_folds,
        outcome_levels  = levels(outcome_f)
    )
}
