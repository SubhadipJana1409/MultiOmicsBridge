#' @title Load and Normalize Microbiome Taxa Table Data
#'
#' @description
#' Imports a microbiome taxa count table and applies either centered
#' log-ratio (CLR) or total sum scaling (TSS) normalization. The result
#' is a \code{SummarizedExperiment} with both raw counts and normalized
#' values stored as named assays.
#'
#' @details
#' Microbiome data is \strong{compositional}: only relative abundances
#' are observed, not absolute counts. Treating compositional data with
#' standard correlation or distance measures leads to spurious results
#' (the Aitchison problem). MultiOmicsBridge applies one of two
#' microbiome-appropriate normalizations:
#'
#' \describe{
#'   \item{CLR (default)}{The centered log-ratio transformation:
#'     \deqn{clr(x_j) = \log(x_j + \delta) - \frac{1}{D}
#'       \sum_{k=1}^{D} \log(x_k + \delta)}
#'     where \eqn{\delta} is the pseudocount and the sum is over all
#'     \eqn{D} taxa. CLR maps compositional data to real space and removes
#'     the unit-sum constraint, enabling Euclidean geometry.}
#'   \item{TSS}{Total sum scaling divides each sample by its library
#'     size, producing relative abundances (proportions). Simpler but
#'     retains the compositional constraint.}
#' }
#'
#' Zero counts are handled by adding a small pseudocount before
#' log-transformation; the default pseudocount of 0.5 is a conservative
#' choice appropriate for sparse 16S data.
#'
#' @param taxa_table A \code{matrix} of non-negative integer counts.
#'   Expected orientation: taxa in \strong{rows}, samples in
#'   \strong{columns} (Bioconductor convention). If samples appear to
#'   be in rows (more rows than columns and row names suggest samples),
#'   the matrix is transposed with a warning.
#' @param col_data An optional \code{data.frame} or \code{DataFrame}
#'   of sample-level metadata. Row names must match the column names of
#'   \code{taxa_table}. Default: \code{NULL}.
#' @param normalization A \code{character(1)}, either \code{"CLR"}
#'   (centered log-ratio, default) or \code{"TSS"} (total sum scaling).
#' @param pseudocount A positive \code{numeric(1)} pseudocount added to
#'   all counts before log-transformation (CLR only). Default: \code{0.5}.
#' @param min_prevalence A \code{numeric(1)} in \code{[0, 1]}. Taxa
#'   present in fewer than this fraction of samples are removed.
#'   Default: \code{0.1} (remove taxa in fewer than 10\% of samples).
#'
#' @return A \code{SummarizedExperiment} with two assays:
#'   \describe{
#'     \item{\code{"counts"}}{Raw integer count matrix (taxa x samples).}
#'     \item{\code{"CLR"} or \code{"TSS"}}{Normalized values.}
#'   }
#'
#' @examples
#' set.seed(42)
#' n_taxa    <- 80
#' n_samples <- 20
#' taxa_table <- matrix(rpois(n_taxa * n_samples, lambda = 30),
#'                      nrow = n_taxa, ncol = n_samples)
#' rownames(taxa_table) <- paste0("Taxon", seq_len(n_taxa))
#' colnames(taxa_table) <- paste0("Sample", seq_len(n_samples))
#'
#' mb_se <- loadMicrobiomeData(taxa_table, normalization = "CLR")
#' mb_se
#' SummarizedExperiment::assayNames(mb_se)
#'
#' @seealso \code{\link{loadHostData}}, \code{\link{matchSamples}},
#'   \code{\link{MultiOmicsBridgeAnalysis}}
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is
#' @export
loadMicrobiomeData <- function(taxa_table,
                               col_data       = NULL,
                               normalization  = c("CLR", "TSS"),
                               pseudocount    = 0.5,
                               min_prevalence = 0.1) {
    normalization <- match.arg(normalization)

    # ── Accept SummarizedExperiment input ────────────────────────────────────
    if (is(taxa_table, "SummarizedExperiment")) {
        if (!"counts" %in% SummarizedExperiment::assayNames(taxa_table))
            stop("'taxa_table' SummarizedExperiment must have a 'counts' assay.")
        if (is.null(col_data))
            col_data <- SummarizedExperiment::colData(taxa_table)
        taxa_table <- SummarizedExperiment::assay(taxa_table, "counts")
    }

    # ── Validation and orientation fix ───────────────────────────────────────
    if (!is.matrix(taxa_table)) taxa_table <- as.matrix(taxa_table)
    if (!is.numeric(taxa_table))
        stop("'taxa_table' must be a numeric matrix of counts.")
    if (any(taxa_table < 0, na.rm = TRUE))
        stop("'taxa_table' must not contain negative values.")

    # Detect potential transposition: if more rows than columns and row names
    # look like sample names (numbered or "Sample_*")
    if (nrow(taxa_table) > ncol(taxa_table) && !is.null(rownames(taxa_table))) {
        sample_like <- grepl("^(Sample|Samp|S)[_\\-]?[0-9]+",
                             rownames(taxa_table)[1L],
                             ignore.case = TRUE)
        if (sample_like) {
            warning("'taxa_table' appears to have samples in rows. ",
                    "Transposing to taxa x samples orientation.")
            taxa_table <- t(taxa_table)
        }
    }

    if (is.null(rownames(taxa_table)))
        stop("'taxa_table' must have taxa names as rownames.")
    if (is.null(colnames(taxa_table)))
        stop("'taxa_table' must have sample names as colnames.")

    if (!is.numeric(pseudocount) || pseudocount <= 0)
        stop("'pseudocount' must be a positive number.")
    if (!is.numeric(min_prevalence) ||
        min_prevalence < 0 || min_prevalence > 1)
        stop("'min_prevalence' must be in [0, 1].")

    # ── Prevalence filtering ─────────────────────────────────────────────────
    n_samples <- ncol(taxa_table)
    prevalence <- rowMeans(taxa_table > 0)
    keep       <- prevalence >= min_prevalence
    if (sum(keep) == 0L)
        stop("All taxa were removed by the 'min_prevalence' filter. ",
             "Try lowering 'min_prevalence'.")
    if (sum(!keep) > 0L)
        message(sprintf(
            "Removing %d low-prevalence taxa (present in < %.0f%% of samples).",
            sum(!keep), min_prevalence * 100))
    taxa_table <- taxa_table[keep, , drop = FALSE]

    message(sprintf("Normalizing %d taxa x %d samples with %s...",
                    nrow(taxa_table), n_samples, normalization))

    # ── Normalization ─────────────────────────────────────────────────────────
    if (normalization == "CLR") {
        normalized <- .clrTransform(taxa_table, pseudocount = pseudocount)
    } else {
        normalized <- .tssTransform(taxa_table)
    }

    # ── Assemble SummarizedExperiment ─────────────────────────────────────────
    cd <- if (is.null(col_data)) {
        S4Vectors::DataFrame(row.names = colnames(taxa_table))
    } else {
        matching <- col_data[match(colnames(taxa_table),
                                   rownames(col_data)), , drop = FALSE]
        rownames(matching) <- colnames(taxa_table)
        S4Vectors::DataFrame(matching)
    }

    assay_list <- list(counts = taxa_table, normalized = normalized)
    names(assay_list)[2L] <- normalization

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = assay_list,
        colData = cd
    )

    message("loadMicrobiomeData complete.")
    se
}
