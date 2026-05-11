#' @title Multi-Omics Biomarker Discovery
#'
#' @description
#' Identifies and ranks host genes and microbial taxa that jointly predict
#' the outcome of interest using two complementary evidence streams: (1)
#' sparse feature loadings from the DIABLO joint dimensionality reduction
#' model, and (2) a cross-omics Spearman correlation network linking host
#' genes to microbial taxa. Features are ranked by their combined loading
#' score and annotated with their maximum cross-omics correlation.
#'
#' @details
#' The biomarker ranking combines:
#' \describe{
#'   \item{Sparse loading score}{The L2 norm of a feature's loadings across
#'     all DIABLO components. Genes/taxa with higher loading scores
#'     contribute more strongly to the latent integration axes.}
#'   \item{Cross-omics correlation}{For each selected host gene, the
#'     maximum absolute Spearman correlation with any selected microbial
#'     taxon (and vice versa). High cross-omics correlation indicates
#'     biologically relevant host-microbe co-variation.}
#' }
#' Hub features — those with both high loading scores and high cross-omics
#' correlations — represent the most credible multi-omics biomarker
#' candidates.
#'
#' @param mae A \code{MultiAssayExperiment} from \code{\link{matchSamples}}.
#' @param dr_result A named \code{list} from \code{\link{jointDimReduction}}.
#' @param n_biomarkers An \code{integer(1)} number of top features to
#'   select per omics layer. Default: \code{50}.
#' @param host_assay A \code{character(1)} assay in host SE for computing
#'   cross-omics correlations. Default: \code{"voom"}.
#' @param mb_assay A \code{character(1)} assay in microbiome SE for
#'   computing cross-omics correlations. Default: \code{"CLR"}.
#'
#' @return A \code{DataFrame} with one row per biomarker and columns:
#'   \describe{
#'     \item{\code{feature}}{Feature name (gene or taxon ID).}
#'     \item{\code{omics_layer}}{Either \code{"host"} or
#'       \code{"microbiome"}.}
#'     \item{\code{loading_score}}{L2 norm of DIABLO loadings across
#'       components.}
#'     \item{\code{rank}}{Within-layer ranking by loading score.}
#'     \item{\code{component}}{DIABLO component with highest absolute
#'       loading.}
#'     \item{\code{max_cross_cor}}{Maximum absolute Spearman correlation
#'       with a feature from the other omics layer.}
#'     \item{\code{top_partner}}{Name of the cross-omics feature with
#'       the highest absolute correlation.}
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
#' bm <- biomarkerDiscovery(mae, dr_res, n_biomarkers = 20)
#' head(as.data.frame(bm))
#'
#' @seealso \code{\link{jointDimReduction}}, \code{\link{plotBiomarkerNetwork}},
#'   \code{\link{MultiOmicsBridgeAnalysis}}
#'
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is
#' @export
biomarkerDiscovery <- function(mae,
                                dr_result,
                                n_biomarkers = 50L,
                                host_assay   = "voom",
                                mb_assay     = "CLR") {
    # ── Validation ────────────────────────────────────────────────────────────
    if (!is(mae, "MultiAssayExperiment"))
        stop("'mae' must be a MultiAssayExperiment.")
    if (!is.list(dr_result) ||
        !all(c("host_loadings", "mb_loadings") %in% names(dr_result)))
        stop("'dr_result' must be the output of jointDimReduction().")

    host_se <- MultiAssayExperiment::experiments(mae)[["host"]]
    mb_se   <- MultiAssayExperiment::experiments(mae)[["microbiome"]]

    host_assay_use <- if (host_assay %in% SummarizedExperiment::assayNames(host_se))
        host_assay else SummarizedExperiment::assayNames(host_se)[1L]
    mb_assay_use   <- if (mb_assay %in% SummarizedExperiment::assayNames(mb_se))
        mb_assay else SummarizedExperiment::assayNames(mb_se)[1L]

    # ── Rank biomarkers from DIABLO loadings ──────────────────────────────────
    message("Ranking biomarkers from DIABLO sparse loadings...")
    bm_df <- .rankBiomarkers(
        host_loadings = dr_result$host_loadings,
        mb_loadings   = dr_result$mb_loadings,
        n_biomarkers  = n_biomarkers
    )

    if (nrow(bm_df) == 0L)
        stop("No biomarkers could be extracted. ",
             "Check that jointDimReduction() ran successfully.")

    # ── Compute cross-omics correlations ──────────────────────────────────────
    message("Computing cross-omics Spearman correlation network...")

    host_mat  <- t(SummarizedExperiment::assay(host_se, host_assay_use))
    mb_mat    <- t(SummarizedExperiment::assay(mb_se,   mb_assay_use))

    # Use only the selected biomarker features
    host_bm_names <- as.character(
        bm_df[["feature"]][bm_df[["omics_layer"]] == "host"])
    mb_bm_names   <- as.character(
        bm_df[["feature"]][bm_df[["omics_layer"]] == "microbiome"])

    host_sel <- host_mat[, intersect(host_bm_names, colnames(host_mat)),
                         drop = FALSE]
    mb_sel   <- mb_mat[,   intersect(mb_bm_names,   colnames(mb_mat)),
                         drop = FALSE]

    cor_mat <- NULL
    if (ncol(host_sel) > 0L && ncol(mb_sel) > 0L && nrow(host_sel) >= 3L) {
        cor_mat <- .crossOmicsCorrelation(host_sel, mb_sel)
    }

    # ── Annotate biomarkers with cross-omics correlations ─────────────────────
    max_cross_cor <- rep(NA_real_, nrow(bm_df))
    top_partner   <- rep(NA_character_, nrow(bm_df))

    if (!is.null(cor_mat)) {
        for (i in seq_len(nrow(bm_df))) {
            feat  <- as.character(bm_df[["feature"]][i])
            layer <- as.character(bm_df[["omics_layer"]][i])

            if (layer == "host" && feat %in% rownames(cor_mat)) {
                row_cors <- abs(cor_mat[feat, ])
                max_idx  <- which.max(row_cors)
                max_cross_cor[i] <- row_cors[max_idx]
                top_partner[i]   <- names(row_cors)[max_idx]
            } else if (layer == "microbiome" && feat %in% colnames(cor_mat)) {
                col_cors <- abs(cor_mat[, feat])
                max_idx  <- which.max(col_cors)
                max_cross_cor[i] <- col_cors[max_idx]
                top_partner[i]   <- names(col_cors)[max_idx]
            }
        }
    }

    bm_df[["max_cross_cor"]] <- max_cross_cor
    bm_df[["top_partner"]]   <- top_partner

    message(sprintf(
        "biomarkerDiscovery: %d host genes, %d microbial taxa selected.",
        sum(bm_df[["omics_layer"]] == "host"),
        sum(bm_df[["omics_layer"]] == "microbiome")
    ))

    bm_df
}
