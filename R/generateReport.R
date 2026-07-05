#' @title Generate a Structured Summary Report of MultiOmicsBridge Results
#'
#' @description
#' Prints a formatted text summary of a \code{\link{MOBResult}} object
#' to the console and optionally saves it to a plain-text file. The
#' report covers all five analysis modules: data dimensions, integration
#' method, top biomarkers, cross-omics correlations, and classifier
#' performance comparison.
#'
#' For a full interactive HTML report, users can render the package
#' vignette template (\code{system.file("vignettes", package =
#' "MultiOmicsBridge")}) with \code{rmarkdown::render()} using their
#' own \code{MOBResult} object.
#'
#' @param result A \code{\link{MOBResult}} object from
#'   \code{\link{MultiOmicsBridgeAnalysis}}.
#' @param file An optional \code{character(1)} file path where the report
#'   should be saved as a plain text file. If \code{NULL} (default), the
#'   report is only printed to the console.
#' @param n_top An \code{integer(1)} number of top biomarkers to list.
#'   Default: \code{10}.
#'
#' @return Invisibly returns a named \code{list} of character vectors,
#'   one per report section.
#'
#' @examples
#' library(S4Vectors)
#' scores <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'     dimnames = list(paste0("S", 1:10), c("Comp1", "Comp2")))
#' bm <- DataFrame(
#'     feature = c("Gene1", "Taxon1", "Gene2"),
#'     omics_layer = c("host", "microbiome", "host"),
#'     loading_score = c(0.9, 0.8, 0.7),
#'     rank = c(1L, 2L, 3L), component = c(1L, 1L, 1L)
#' )
#' cr <- list(
#'     host_only       = list(mean_auc = 0.82, sd_auc = 0.05,
#'                            fold_auc = c(0.78, 0.84, 0.83)),
#'     microbiome_only = list(mean_auc = 0.75, sd_auc = 0.06,
#'                            fold_auc = c(0.70, 0.79, 0.76)),
#'     joint           = list(mean_auc = 0.94, sd_auc = 0.03,
#'                            fold_auc = c(0.92, 0.95, 0.94)),
#'     cv_folds        = 3L,
#'     outcome_levels  = c("ctrl", "treat")
#' )
#' obj <- MOBResult(scores, list(), bm, cr,
#'     params = list(integration_method = "DIABLO",
#'                   n_components = 2L,
#'                   outcome_levels = c("ctrl", "treat"),
#'                   cv_folds = 3L))
#' generateReport(obj, n_top = 3)
#'
#' @seealso \code{\link{MultiOmicsBridgeAnalysis}}, \code{\link{MOBResult}}
#'
#' @importFrom methods is
#' @importFrom utils head
#' @export
generateReport <- function(result, file = NULL, n_top = 10L) {
    if (!is(result, "MOBResult"))
        stop("'result' must be a MOBResult object.")

    p    <- result@params
    bm   <- as.data.frame(biomarkers(result))
    cr   <- performance(result)
    sc   <- integrationScores(result)

    lines <- list()

    sep <- paste(rep("=", 60), collapse = "")
    dash <- paste(rep("-", 60), collapse = "")

    # ── Header ────────────────────────────────────────────────────────────────
    lines$header <- c(
        sep,
        "  MultiOmicsBridge Analysis Report",
        sprintf("  Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
        sep
    )

    # ── Module 1: Data Summary ────────────────────────────────────────────────
    lines$data_summary <- c(
        "",
        "MODULE 1: Data Summary",
        dash,
        sprintf("  Paired samples      : %d", nrow(sc)),
        sprintf("  Latent components   : %d", ncol(sc)),
        sprintf("  Integration method  : %s",
                p$integration_method %||% "DIABLO"),
        sprintf("  Outcome             : %s",
                paste(p$outcome_levels %||% "unknown", collapse = " vs ")),
        sprintf("  Host features used  : %s per component",
                as.character(p$n_features_host %||% "N/A")),
        sprintf("  MB features used    : %s per component",
                as.character(p$n_features_mb %||% "N/A"))
    )

    # ── Module 3: Biomarker Summary ───────────────────────────────────────────
    n_host  <- if (nrow(bm) > 0) sum(bm$omics_layer == "host") else 0L
    n_mb    <- if (nrow(bm) > 0) sum(bm$omics_layer == "microbiome") else 0L

    top_bm <- utils::head(bm[order(bm$loading_score, decreasing = TRUE), ], n_top)

    bm_lines <- c(
        "",
        "MODULE 3: Multi-Omics Biomarkers",
        dash,
        sprintf("  Total biomarkers    : %d (%d host, %d microbiome)",
                nrow(bm), n_host, n_mb)
    )

    if (nrow(top_bm) > 0L) {
        bm_lines <- c(bm_lines,
                      sprintf("  Top %d biomarkers (by loading score):", n_top),
                      sprintf("  %4s %-25s %-12s %s",
                              "Rank", "Feature", "Layer", "Loading"))
        for (i in seq_len(nrow(top_bm))) {
            bm_lines <- c(bm_lines,
                          sprintf("  %4d %-25s %-12s %.4f",
                                  top_bm$rank[i],
                                  top_bm$feature[i],
                                  top_bm$omics_layer[i],
                                  top_bm$loading_score[i]))
        }
        if (!is.null(bm$max_cross_cor)) {
            hub <- top_bm[!is.na(top_bm$max_cross_cor), ]
            if (nrow(hub) > 0L) {
                hub <- hub[order(hub$max_cross_cor, decreasing = TRUE), ]
                bm_lines <- c(bm_lines, "",
                              "  Top cross-omics hubs (highest partner correlation):")
                for (i in seq_len(min(5L, nrow(hub)))) {
                    bm_lines <- c(bm_lines,
                                  sprintf("  %-25s -> %-25s (r=%.3f)",
                                          hub$feature[i],
                                          hub$top_partner[i] %||% "NA",
                                          hub$max_cross_cor[i]))
                }
            }
        }
    }
    lines$biomarkers <- bm_lines

    # ── Module 4: Classifier Performance ─────────────────────────────────────
    clf_lines <- c(
        "",
        "MODULE 4: Diagnostic Classifier Performance",
        dash,
        sprintf("  Cross-validation    : %d-fold",
                cr$cv_folds %||% p$cv_folds %||% 5L)
    )

    model_labels <- c(host_only       = "Host RNA-seq only",
                      microbiome_only = "Microbiome only",
                      joint           = "Joint (multi-omics)")

    best_auc  <- -Inf
    best_model <- "joint"

    for (nm in c("host_only", "microbiome_only", "joint")) {
        r <- cr[[nm]]
        if (is.null(r)) next
        clf_lines <- c(clf_lines,
                       sprintf("  %-22s  AUC = %.3f +/- %.3f",
                               model_labels[nm],
                               r$mean_auc %||% NA,
                               r$sd_auc   %||% 0))
        if (!is.na(r$mean_auc) && r$mean_auc > best_auc) {
            best_auc   <- r$mean_auc
            best_model <- nm
        }
    }

    if (!is.null(cr$joint) && !is.null(cr$host_only)) {
        delta <- (cr$joint$mean_auc %||% 0) - (cr$host_only$mean_auc %||% 0)
        clf_lines <- c(clf_lines, "",
                       sprintf("  Multi-omics gain vs host-only: +%.3f AUC", delta))
    }
    lines$classifier <- clf_lines

    # ── Footer ────────────────────────────────────────────────────────────────
    lines$footer <- c("", sep, "  End of report", sep)

    # ── Print ─────────────────────────────────────────────────────────────────
    all_lines <- unlist(lines)
    message(paste(all_lines, collapse = "\n"))

    # ── Optional file save ────────────────────────────────────────────────────
    if (!is.null(file)) {
        writeLines(all_lines, con = file)
        message("Report saved to: ", file)
    }

    invisible(lines)
}
