#' @title Joint Dimensionality Reduction via Sparse Multi-Block PLS-DA
#'
#' @description
#' Performs joint dimensionality reduction across host transcriptomics
#' and gut microbiome data using DIABLO (Data Integration Analysis for
#' Biomarker discovery using Latent cOmponents), a sparse multi-block
#' PLS-DA implemented in \code{mixOmics}. DIABLO simultaneously identifies
#' correlated features across data blocks while discriminating between
#' outcome groups, enforcing sparsity so only the most informative
#' features contribute to each latent component.
#'
#' @details
#' The DIABLO model is fitted as:
#' \deqn{\max \, \text{Cov}(t_\text{host}, t_\text{mb})}
#' subject to \eqn{||w_\text{host}||_1 \leq \lambda_\text{host}} and
#' \eqn{||w_\text{mb}||_1 \leq \lambda_\text{mb}}, where
#' \eqn{t} are the sample scores (variates) and \eqn{w} are the
#' sparse feature weights (loadings).
#'
#' The number of features retained per component is controlled by
#' \code{n_features_host} and \code{n_features_mb}. By default, a
#' design matrix connecting all blocks with moderate correlation
#' (\code{design_off_diag = 0.1}) is used, which prioritizes outcome
#' discrimination over cross-block correlation.
#'
#' @param mae A \code{MultiAssayExperiment} from \code{\link{matchSamples}}
#'   containing \code{"host"} and \code{"microbiome"} experiments.
#' @param outcome A \code{character} or \code{factor} vector of outcome
#'   labels, one per sample (in the same order as \code{colnames(mae)}).
#'   Must have exactly 2 or more levels.
#' @param n_components An \code{integer(1)} number of latent components
#'   to extract. Default: \code{2}.
#' @param n_features_host An \code{integer(1)} number of host genes to
#'   retain per component (sparse selection). Default: \code{50}.
#' @param n_features_mb An \code{integer(1)} number of microbial taxa to
#'   retain per component. Default: \code{20}.
#' @param design_off_diag A \code{numeric(1)} in \code{[0, 1]} controlling
#'   the off-diagonal entries of the DIABLO design matrix. Higher values
#'   emphasize cross-block correlation; lower values prioritize outcome
#'   discrimination. Default: \code{0.1}.
#' @param host_assay A \code{character(1)} name of the assay in the host
#'   \code{SummarizedExperiment} to use as input. Default: \code{"voom"}.
#' @param mb_assay A \code{character(1)} name of the assay in the
#'   microbiome \code{SummarizedExperiment} to use as input. Default:
#'   \code{"CLR"}.
#' @param min_variance A \code{numeric(1)} minimum variance (as a fraction
#'   of total row variance) for genes/taxa to be retained before DIABLO.
#'   Default: \code{0} (no filtering).
#'
#' @return A named \code{list} with elements:
#'   \describe{
#'     \item{\code{scores}}{A \code{matrix} (samples x components) of
#'       integrated sample scores (from the host block variate).}
#'     \item{\code{host_loadings}}{A \code{matrix} (genes x components)
#'       of sparse host feature loadings.}
#'     \item{\code{mb_loadings}}{A \code{matrix} (taxa x components)
#'       of sparse microbiome feature loadings.}
#'     \item{\code{explained_variance}}{A named \code{numeric} vector of
#'       explained variance per component.}
#'     \item{\code{diablo_object}}{The full DIABLO result object from
#'       \code{mixOmics::block.splsda}.}
#'     \item{\code{outcome}}{The outcome factor used.}
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
#'
#' outcome <- rep(c("ctrl", "treat"), each = 10)
#' dr_res  <- jointDimReduction(mae, outcome = outcome,
#'                               n_components = 2,
#'                               n_features_host = 30,
#'                               n_features_mb   = 15)
#' dim(dr_res$scores)
#' head(dr_res$explained_variance)
#'
#' @seealso \code{\link{biomarkerDiscovery}}, \code{\link{plotIntegration}},
#'   \code{\link{MultiOmicsBridgeAnalysis}}
#'
#' @importFrom MultiAssayExperiment experiments colData
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom mixOmics block.splsda
#' @importFrom methods is
#' @importFrom stats quantile
#' @export
jointDimReduction <- function(mae,
                               outcome,
                               n_components     = 2L,
                               n_features_host  = 50L,
                               n_features_mb    = 20L,
                               design_off_diag  = 0.1,
                               host_assay       = "voom",
                               mb_assay         = "CLR",
                               min_variance     = 0) {
    # ── Validation ────────────────────────────────────────────────────────────
    if (!is(mae, "MultiAssayExperiment"))
        stop("'mae' must be a MultiAssayExperiment object from matchSamples().")

    exp_names <- names(MultiAssayExperiment::experiments(mae))
    if (!all(c("host", "microbiome") %in% exp_names))
        stop("'mae' must contain 'host' and 'microbiome' experiments. ",
             "Found: ", paste(exp_names, collapse = ", "))

    n_samples <- nrow(MultiAssayExperiment::colData(mae))
    if (length(outcome) != n_samples)
        stop("Length of 'outcome' (", length(outcome), ") must equal ",
             "number of samples (", n_samples, ").")

    outcome_f <- factor(outcome)
    if (nlevels(outcome_f) < 2L)
        stop("'outcome' must have at least 2 levels.")
    if (!is.numeric(n_components) || n_components < 1L)
        stop("'n_components' must be a positive integer.")
    if (!is.numeric(design_off_diag) ||
        design_off_diag < 0 || design_off_diag > 1)
        stop("'design_off_diag' must be in [0, 1].")

    # ── Extract assay matrices ────────────────────────────────────────────────
    host_se <- MultiAssayExperiment::experiments(mae)[["host"]]
    mb_se   <- MultiAssayExperiment::experiments(mae)[["microbiome"]]

    host_assay_use <- if (host_assay %in% SummarizedExperiment::assayNames(host_se))
        host_assay else {
            warning("Assay '", host_assay, "' not found in host SE. ",
                    "Using first available assay.")
            SummarizedExperiment::assayNames(host_se)[1L]
        }

    mb_assay_use <- if (mb_assay %in% SummarizedExperiment::assayNames(mb_se))
        mb_assay else {
            warning("Assay '", mb_assay, "' not found in microbiome SE. ",
                    "Using first available assay.")
            SummarizedExperiment::assayNames(mb_se)[1L]
        }

    # DIABLO expects: samples x features
    host_mat <- t(SummarizedExperiment::assay(host_se, host_assay_use))
    mb_mat   <- t(SummarizedExperiment::assay(mb_se,   mb_assay_use))

    # ── Optional variance filtering ───────────────────────────────────────────
    if (min_variance > 0) {
        host_var <- apply(host_mat, 2L, var)
        mb_var   <- apply(mb_mat,   2L, var)
        host_mat <- host_mat[, host_var >= quantile(host_var, min_variance),
                             drop = FALSE]
        mb_mat   <- mb_mat[, mb_var >= quantile(mb_var, min_variance),
                           drop = FALSE]
        message(sprintf("After variance filter: %d host genes, %d taxa.",
                        ncol(host_mat), ncol(mb_mat)))
    }

    n_features_host <- min(n_features_host, ncol(host_mat))
    n_features_mb   <- min(n_features_mb,   ncol(mb_mat))

    message(sprintf(
        "Running DIABLO: %d samples | %d host features | %d mb features | %d components",
        nrow(host_mat), ncol(host_mat), ncol(mb_mat), n_components
    ))

    # ── Build DIABLO design matrix ────────────────────────────────────────────
    # 2x2 design: rows/cols = blocks (host, microbiome)
    # 0 on diagonal, design_off_diag off-diagonal
    design <- matrix(design_off_diag, nrow = 2L, ncol = 2L,
                     dimnames = list(c("host","microbiome"),
                                     c("host","microbiome")))
    diag(design) <- 0

    # ── Run DIABLO ────────────────────────────────────────────────────────────
    diablo_res <- tryCatch(
        mixOmics::block.splsda(
            X      = list(host       = host_mat,
                          microbiome = mb_mat),
            Y      = outcome_f,
            ncomp  = n_components,
            keepX  = list(
                host       = rep(n_features_host, n_components),
                microbiome = rep(n_features_mb,   n_components)
            ),
            design = design
        ),
        error = function(e) stop("DIABLO failed: ", conditionMessage(e))
    )

    # ── Extract scores and loadings ───────────────────────────────────────────
    scores        <- .extractScores(diablo_res, "host")
    host_loadings <- .extractLoadings(diablo_res, "host")
    mb_loadings   <- .extractLoadings(diablo_res, "microbiome")

    # Explained variance: proportion of outcome variance explained
    exp_var <- tryCatch({
        ev <- diablo_res[["prop_expl_var"]]
        if (!is.null(ev[["host"]])) ev[["host"]] else rep(NA_real_, n_components)
    }, error = function(e) rep(NA_real_, n_components))
    names(exp_var) <- paste0("Comp", seq_len(n_components))

    message("jointDimReduction complete.")
    list(
        scores            = scores,
        host_loadings     = host_loadings,
        mb_loadings       = mb_loadings,
        explained_variance = exp_var,
        diablo_object     = diablo_res,
        outcome           = outcome_f
    )
}
