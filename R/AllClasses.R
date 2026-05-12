#' @title MOBResult: Multi-Omics Bridge Result Container
#'
#' @description
#' An S4 class storing the output of \code{\link{MultiOmicsBridgeAnalysis}}.
#' Slots hold integrated sample scores from joint dimensionality reduction,
#' per-layer feature loadings, a ranked multi-omics biomarker table,
#' diagnostic classifier performance metrics, and analysis parameters.
#'
#' @slot integratedScores A \code{matrix} of integrated sample scores with
#'   rows = samples and columns = latent components. Produced by DIABLO.
#' @slot featureLoadings A named \code{list} with two elements:
#'   \code{host} (genes x components) and \code{microbiome} (taxa x
#'   components) containing the sparse feature loading matrices.
#' @slot biomarkerTable A \code{DataFrame} with one row per selected
#'   biomarker containing columns \code{feature}, \code{omics_layer},
#'   \code{loading_score}, \code{rank}, and \code{component}.
#' @slot classifierResults A named \code{list} with elements
#'   \code{host_only}, \code{microbiome_only}, and \code{joint}, each
#'   containing cross-validated AUC-ROC and fold-level performance metrics.
#' @slot params A \code{list} of analysis parameters including
#'   \code{integration_method}, \code{n_components}, \code{n_biomarkers},
#'   \code{cv_folds}, and \code{outcome_levels}.
#'
#' @exportClass MOBResult
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new is
setClass(
    "MOBResult",
    representation(
        integratedScores  = "matrix",
        featureLoadings   = "list",
        biomarkerTable    = "DataFrame",
        classifierResults = "list",
        params            = "list"
    )
)

# ── Constructor ───────────────────────────────────────────────────────────────

#' @title Constructor for MOBResult
#'
#' @description Create a new \code{MOBResult} object.
#'
#' @param integratedScores A \code{matrix} of integrated sample scores
#'   (samples x components).
#' @param featureLoadings A named \code{list} with \code{host} and
#'   \code{microbiome} loading matrices.
#' @param biomarkerTable A \code{DataFrame} of ranked multi-omics
#'   biomarkers.
#' @param classifierResults A named \code{list} of classifier performance
#'   metrics.
#' @param params A \code{list} of analysis parameters.
#'
#' @return A \code{MOBResult} object.
#'
#' @examples
#' library(S4Vectors)
#' scores <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'                  dimnames = list(paste0("S", 1:10), c("Comp1","Comp2")))
#' loadings <- list(
#'     host       = matrix(rnorm(10), nrow = 5, ncol = 2,
#'                         dimnames = list(paste0("G", 1:5), c("Comp1","Comp2"))),
#'     microbiome = matrix(rnorm(6), nrow = 3, ncol = 2,
#'                         dimnames = list(paste0("T", 1:3), c("Comp1","Comp2")))
#' )
#' bm <- DataFrame(
#'     feature      = c("G1", "T1"),
#'     omics_layer  = c("host", "microbiome"),
#'     loading_score = c(0.8, 0.6),
#'     rank         = c(1L, 2L),
#'     component    = c(1L, 1L)
#' )
#' obj <- MOBResult(
#'     integratedScores  = scores,
#'     featureLoadings   = loadings,
#'     biomarkerTable    = bm,
#'     classifierResults = list(),
#'     params = list(integration_method = "DIABLO")
#' )
#' obj
#'
#' @export
MOBResult <- function(integratedScores, featureLoadings, biomarkerTable,
                      classifierResults, params = list()) {
    new("MOBResult",
        integratedScores  = integratedScores,
        featureLoadings   = featureLoadings,
        biomarkerTable    = biomarkerTable,
        classifierResults = classifierResults,
        params            = params)
}

# ── Generics ──────────────────────────────────────────────────────────────────

#' @title Accessor for biomarker table in a MOBResult
#'
#' @description Returns the ranked multi-omics biomarker \code{DataFrame}
#'   from a \code{MOBResult} object.
#'
#' @param x A \code{MOBResult} object.
#' @param ... Additional arguments (not used).
#'
#' @return A \code{DataFrame} with columns \code{feature},
#'   \code{omics_layer}, \code{loading_score}, \code{rank}, and
#'   \code{component}.
#'
#' @examples
#' library(S4Vectors)
#' bm <- DataFrame(feature = c("G1","T1"), omics_layer = c("host","microbiome"),
#'                 loading_score = c(0.8,0.6), rank = c(1L,2L), component = c(1L,1L))
#' obj <- MOBResult(matrix(rnorm(20), 10, 2), list(), bm, list())
#' biomarkers(obj)
#'
#' @export
setGeneric("biomarkers", function(x, ...) standardGeneric("biomarkers"))

#' @title Accessor for classifier performance in a MOBResult
#'
#' @description Returns the cross-validated classifier performance metrics
#'   from a \code{MOBResult} object.
#'
#' @param x A \code{MOBResult} object.
#' @param ... Additional arguments (not used).
#'
#' @return A named \code{list} with elements \code{host_only},
#'   \code{microbiome_only}, and \code{joint}, each containing
#'   \code{mean_auc}, \code{sd_auc}, and \code{fold_auc}.
#'
#' @examples
#' library(S4Vectors)
#' cr <- list(host_only = list(mean_auc = 0.85),
#'            microbiome_only = list(mean_auc = 0.78),
#'            joint = list(mean_auc = 0.92))
#' obj <- MOBResult(matrix(rnorm(20), 10, 2), list(),
#'                  DataFrame(), cr)
#' performance(obj)
#'
#' @export
setGeneric("performance", function(x, ...) standardGeneric("performance"))

#' @title Accessor for integrated sample scores in a MOBResult
#'
#' @description Returns the matrix of integrated sample scores (samples x
#'   latent components) from a \code{MOBResult} object.
#'
#' @param x A \code{MOBResult} object.
#' @param ... Additional arguments (not used).
#'
#' @return A \code{matrix} with rows = samples and columns = latent
#'   components.
#'
#' @examples
#' library(S4Vectors)
#' scores <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'                  dimnames = list(paste0("S",1:10), c("Comp1","Comp2")))
#' obj <- MOBResult(scores, list(), DataFrame(), list())
#' integrationScores(obj)
#'
#' @export
setGeneric("integrationScores",
           function(x, ...) standardGeneric("integrationScores"))

#' @title Accessor for feature loadings in a MOBResult
#'
#' @description Returns the named list of per-layer feature loading
#'   matrices from a \code{MOBResult} object.
#'
#' @param x A \code{MOBResult} object.
#' @param ... Additional arguments (not used).
#'
#' @return A named \code{list} with elements \code{host} and
#'   \code{microbiome}, each a matrix of genes/taxa x components.
#'
#' @examples
#' library(S4Vectors)
#' fl <- list(host = matrix(rnorm(10), 5, 2), microbiome = matrix(rnorm(6), 3, 2))
#' obj <- MOBResult(matrix(rnorm(20), 10, 2), fl, DataFrame(), list())
#' featureLoadings(obj)
#'
#' @export
setGeneric("featureLoadings",
           function(x, ...) standardGeneric("featureLoadings"))

# ── Methods ───────────────────────────────────────────────────────────────────

#' @title Show method for MOBResult
#'
#' @description Prints a compact summary of a \code{MOBResult} object,
#'   including the top biomarkers, integration method, and classifier
#'   AUC-ROC values.
#'
#' @param object A \code{MOBResult} object.
#'
#' @return Invisibly returns \code{object}.
#'
#' @importFrom methods show
#' @export
setMethod("show", "MOBResult", function(object) {
    p    <- object@params
    bm   <- object@biomarkerTable
    cr   <- object@classifierResults
    sc   <- object@integratedScores

    n_host <- if ("omics_layer" %in% names(bm))
        sum(bm[["omics_layer"]] == "host") else NA_integer_
    n_mb   <- if ("omics_layer" %in% names(bm))
        sum(bm[["omics_layer"]] == "microbiome") else NA_integer_

    cat("MOBResult\n")
    cat("  Integration    :",
        if (!is.null(p$integration_method)) p$integration_method
        else "not set", "\n")
    cat("  Samples        :", nrow(sc), "\n")
    cat("  Components     :", ncol(sc), "\n")
    cat("  Biomarkers     :", nrow(bm),
        sprintf("(%d host, %d microbiome)", n_host, n_mb), "\n")
    cat("  -- Classifier AUC (mean +/- SD) -----------\n")
    for (nm in c("host_only", "microbiome_only", "joint")) {
        r <- cr[[nm]]
        if (!is.null(r)) {
            cat(sprintf("  %-16s %.3f +/- %.3f\n",
                        paste0(nm, ":"),
                        r$mean_auc %||% NA,
                        r$sd_auc   %||% NA))
        }
    }
    if (!is.null(p$outcome_levels))
        cat("  Outcome levels :",
            paste(p$outcome_levels, collapse = " vs "), "\n")
    invisible(object)
})

#' @rdname biomarkers
#' @export
setMethod("biomarkers", "MOBResult", function(x, ...) x@biomarkerTable)

#' @rdname performance
#' @export
setMethod("performance", "MOBResult", function(x, ...) x@classifierResults)

#' @rdname integrationScores
#' @export
setMethod("integrationScores", "MOBResult",
          function(x, ...) x@integratedScores)

#' @rdname featureLoadings
#' @export
setMethod("featureLoadings", "MOBResult",
          function(x, ...) x@featureLoadings)

# ── Utility ───────────────────────────────────────────────────────────────────

#' Null-coalescing operator (internal)
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b
