#' @title Load and Normalize Host Bulk RNA-seq Data
#'
#' @description
#' Imports a bulk RNA-seq count matrix and applies TMM normalization
#' followed by limma-voom precision weighting. The result is a
#' \code{SummarizedExperiment} with both the raw counts and
#' voom-transformed log2-CPM values stored as named assays.
#'
#' @details
#' The normalization pipeline:
#' \enumerate{
#'   \item Construct a \code{DGEList} from raw counts.
#'   \item Apply TMM (trimmed mean of M-values) normalization via
#'     \code{edgeR::calcNormFactors} to remove compositional bias
#'     between libraries.
#'   \item Apply \code{limma::voom} to compute log2-CPM values with
#'     precision weights that model the mean-variance trend. These
#'     weights are used downstream by \code{\link{diagnosticClassifier}}
#'     and \code{\link{jointDimReduction}}.
#' }
#'
#' Genes with very low expression are not filtered here; gene filtering
#' is performed within \code{\link{jointDimReduction}} and
#' \code{\link{diagnosticClassifier}} based on expression in the matched
#' multi-omics dataset. Users may pre-filter genes if desired.
#'
#' @param counts A \code{matrix} of non-negative integer counts with
#'   genes in rows and samples in columns. Alternatively, a
#'   \code{SummarizedExperiment} whose \code{"counts"} assay is used.
#' @param col_data An optional \code{data.frame} or \code{DataFrame}
#'   of sample-level metadata. Row names must match the column names of
#'   \code{counts}. If \code{NULL}, an empty \code{DataFrame} is used.
#'   Default: \code{NULL}.
#' @param min_count A \code{numeric(1)} minimum total count per gene.
#'   Genes with row sums below this threshold are removed before
#'   normalization. Default: \code{1} (removes all-zero genes only).
#' @param assay_name A \code{character(1)} name for the voom-normalized
#'   assay in the output \code{SummarizedExperiment}. Default:
#'   \code{"voom"}.
#'
#' @return A \code{SummarizedExperiment} with two assays:
#'   \describe{
#'     \item{\code{"counts"}}{Raw integer count matrix.}
#'     \item{\code{"voom"}}{Voom-transformed log2-CPM matrix with
#'       TMM library size normalization applied.}
#'   }
#'   Sample metadata (if provided) is stored in \code{colData}.
#'
#' @examples
#' set.seed(42)
#' n_genes   <- 200
#' n_samples <- 20
#' counts <- matrix(rpois(n_genes * n_samples, lambda = 150),
#'                  nrow = n_genes, ncol = n_samples)
#' rownames(counts) <- paste0("Gene", seq_len(n_genes))
#' colnames(counts) <- paste0("Sample", seq_len(n_samples))
#'
#' col_data <- data.frame(
#'     condition = rep(c("ctrl", "treat"), each = 10),
#'     row.names = colnames(counts)
#' )
#'
#' host_se <- loadHostData(counts, col_data = col_data)
#' host_se
#' SummarizedExperiment::assayNames(host_se)
#'
#' @seealso \code{\link{loadMicrobiomeData}}, \code{\link{matchSamples}},
#'   \code{\link{MultiOmicsBridgeAnalysis}}
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames
#' @importFrom S4Vectors DataFrame
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom limma voom
#' @importFrom methods is
#' @export
loadHostData <- function(counts,
                         col_data  = NULL,
                         min_count = 1,
                         assay_name = "voom") {
    # ── Accept SummarizedExperiment input ────────────────────────────────────
    if (is(counts, "SummarizedExperiment")) {
        if (!"counts" %in% SummarizedExperiment::assayNames(counts))
            stop("'counts' SummarizedExperiment must have a 'counts' assay.")
        if (is.null(col_data))
            col_data <- SummarizedExperiment::colData(counts)
        counts <- SummarizedExperiment::assay(counts, "counts")
    }

    # ── Validation ────────────────────────────────────────────────────────────
    if (!is.matrix(counts)) {
        counts <- as.matrix(counts)
        if (!is.numeric(counts))
            stop("'counts' must be a numeric matrix of raw counts.")
    }
    if (is.null(rownames(counts)))
        stop("'counts' must have gene names as rownames.")
    if (is.null(colnames(counts)))
        stop("'counts' must have sample names as colnames.")
    if (any(counts < 0, na.rm = TRUE))
        stop("'counts' must not contain negative values.")
    if (!is.numeric(min_count) || min_count < 0)
        stop("'min_count' must be a non-negative number.")

    # ── Filter lowly expressed genes ─────────────────────────────────────────
    row_sums <- rowSums(counts)
    keep     <- row_sums >= min_count
    if (sum(keep) == 0L)
        stop("All genes were removed by the 'min_count' filter.")
    if (sum(!keep) > 0L)
        message(sprintf("Removing %d zero/low-count genes (min_count = %g).",
                        sum(!keep), min_count))
    counts <- counts[keep, , drop = FALSE]

    # ── TMM normalization + voom ─────────────────────────────────────────────
    message(sprintf("Normalizing %d genes x %d samples with TMM + voom...",
                    nrow(counts), ncol(counts)))
    voom_res <- .voomNormalize(counts, design = NULL)

    # ── Assemble SummarizedExperiment ─────────────────────────────────────────
    cd <- if (is.null(col_data)) {
        S4Vectors::DataFrame(row.names = colnames(counts))
    } else {
        # Align col_data rownames to sample order
        if (!all(colnames(counts) %in% rownames(col_data)))
            warning("Not all sample names in 'counts' are in 'col_data'. ",
                    "Proceeding with empty metadata for missing samples.")
        matching <- col_data[match(colnames(counts), rownames(col_data)), ,
                             drop = FALSE]
        rownames(matching) <- colnames(counts)
        S4Vectors::DataFrame(matching)
    }

    assay_list <- list(
        counts = counts,
        voom   = voom_res$E
    )
    names(assay_list)[2L] <- assay_name

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = assay_list,
        colData = cd
    )

    message("loadHostData complete.")
    se
}
