#' @title Match Paired Samples Across Host and Microbiome Datasets
#'
#' @description
#' Identifies samples present in both host and microbiome
#' \code{SummarizedExperiment} objects, subsets both to the common
#' samples, and assembles a \code{MultiAssayExperiment} (MAE) that
#' serves as the primary input for downstream analysis functions.
#'
#' @details
#' In paired multi-omics studies, not all samples necessarily have both
#' data types due to sequencing failures, QC exclusions, or study design.
#' \code{matchSamples} transparently reports how many samples are retained
#' and warns if fewer than \code{min_paired} paired samples are found.
#'
#' The output \code{MultiAssayExperiment} stores the host and microbiome
#' \code{SummarizedExperiment} objects under the names \code{"host"} and
#' \code{"microbiome"} respectively, with a unified \code{colData} drawn
#' from the host sample metadata.
#'
#' @param host_se A \code{SummarizedExperiment} from
#'   \code{\link{loadHostData}}, with samples as columns.
#' @param mb_se A \code{SummarizedExperiment} from
#'   \code{\link{loadMicrobiomeData}}, with samples as columns.
#' @param sample_col_host A \code{character(1)} naming the column in
#'   \code{colData(host_se)} that contains the sample identifier. If
#'   \code{NULL} (default), uses \code{colnames(host_se)}.
#' @param sample_col_mb A \code{character(1)} naming the column in
#'   \code{colData(mb_se)} that contains the sample identifier. If
#'   \code{NULL} (default), uses \code{colnames(mb_se)}.
#' @param min_paired An \code{integer(1)} minimum number of paired
#'   samples required to proceed. Default: \code{5}.
#'
#' @return A \code{MultiAssayExperiment} with two experiments:
#'   \describe{
#'     \item{\code{"host"}}{Subset of \code{host_se} for paired samples.}
#'     \item{\code{"microbiome"}}{Subset of \code{mb_se} for paired
#'       samples.}
#'   }
#'   The \code{colData} of the MAE is taken from \code{host_se} for the
#'   paired samples.
#'
#' @examples
#' set.seed(42)
#' # Host data: 200 genes, 20 samples
#' host_counts <- matrix(rpois(200 * 20, 150),
#'                       nrow = 200, ncol = 20,
#'                       dimnames = list(paste0("Gene", 1:200),
#'                                       paste0("Sample", 1:20)))
#' host_se <- loadHostData(host_counts)
#'
#' # Microbiome data: 50 taxa, 18 samples (2 missing)
#' mb_counts <- matrix(rpois(50 * 18, 30),
#'                     nrow = 50, ncol = 18,
#'                     dimnames = list(paste0("Taxon", 1:50),
#'                                     paste0("Sample", 1:18)))
#' mb_se <- loadMicrobiomeData(mb_counts)
#'
#' mae <- matchSamples(host_se, mb_se, min_paired = 5)
#' mae
#'
#' @seealso \code{\link{loadHostData}}, \code{\link{loadMicrobiomeData}},
#'   \code{\link{MultiOmicsBridgeAnalysis}}
#'
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is
#' @export
matchSamples <- function(host_se,
                         mb_se,
                         sample_col_host = NULL,
                         sample_col_mb   = NULL,
                         min_paired      = 5L) {
    # ── Input validation ──────────────────────────────────────────────────────
    if (!is(host_se, "SummarizedExperiment"))
        stop("'host_se' must be a SummarizedExperiment object.")
    if (!is(mb_se, "SummarizedExperiment"))
        stop("'mb_se' must be a SummarizedExperiment object.")
    if (!is.numeric(min_paired) || min_paired < 2L)
        stop("'min_paired' must be at least 2.")

    # ── Extract sample IDs ────────────────────────────────────────────────────
    host_ids <- if (!is.null(sample_col_host) &&
                    sample_col_host %in% names(SummarizedExperiment::colData(host_se))) {
        as.character(SummarizedExperiment::colData(host_se)[[sample_col_host]])
    } else {
        colnames(host_se)
    }

    mb_ids <- if (!is.null(sample_col_mb) &&
                  sample_col_mb %in% names(SummarizedExperiment::colData(mb_se))) {
        as.character(SummarizedExperiment::colData(mb_se)[[sample_col_mb]])
    } else {
        colnames(mb_se)
    }

    # ── Find common samples ───────────────────────────────────────────────────
    common_ids <- intersect(host_ids, mb_ids)
    n_common   <- length(common_ids)

    message(sprintf(
        "matchSamples: %d host | %d microbiome | %d paired samples retained.",
        length(host_ids), length(mb_ids), n_common
    ))

    if (n_common == 0L)
        stop("No common sample names found between 'host_se' and 'mb_se'. ",
             "Ensure column names (or the specified 'sample_col' columns) ",
             "match between the two objects.")

    if (n_common < min_paired)
        stop(sprintf(
            "Only %d paired samples found. At least %d required. ",
            n_common, min_paired),
            "Check sample name consistency or lower 'min_paired'.")

    if (n_common < length(host_ids) || n_common < length(mb_ids)) {
        host_only <- setdiff(host_ids, mb_ids)
        mb_only   <- setdiff(mb_ids, host_ids)
        if (length(host_only) > 0L)
            message(sprintf("  %d host-only samples excluded: %s",
                            length(host_only),
                            paste(utils::head(host_only, 5L), collapse = ", ")))
        if (length(mb_only) > 0L)
            message(sprintf("  %d microbiome-only samples excluded: %s",
                            length(mb_only),
                            paste(utils::head(mb_only, 5L), collapse = ", ")))
    }

    # ── Subset both SEs to common samples ────────────────────────────────────
    host_idx <- match(common_ids, host_ids)
    mb_idx   <- match(common_ids, mb_ids)

    host_sub <- host_se[, host_idx]
    mb_sub   <- mb_se[,   mb_idx]

    # Standardize colnames to common_ids
    colnames(host_sub) <- common_ids
    colnames(mb_sub)   <- common_ids

    # ── Build MultiAssayExperiment ────────────────────────────────────────────
    primary_cd <- SummarizedExperiment::colData(host_sub)

    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(host       = host_sub,
                           microbiome = mb_sub),
        colData     = primary_cd
    )

    mae
}
