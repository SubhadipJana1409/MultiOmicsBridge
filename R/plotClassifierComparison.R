#' @title ROC Curve Comparison of Host-Only, Microbiome-Only, and Joint Classifiers
#'
#' @description
#' Produces overlaid ROC (Receiver Operating Characteristic) curves
#' comparing the three diagnostic classifier configurations — host
#' transcriptomics only, gut microbiome only, and the joint multi-omics
#' model — on the same axes. The AUC-ROC values are annotated on the
#' plot, making the multi-omics advantage immediately visible. A
#' bar chart of mean cross-validated AUC values is also available via
#' \code{type = "bar"}.
#'
#' @param result A \code{\link{MOBResult}} object from
#'   \code{\link{MultiOmicsBridgeAnalysis}}.
#' @param type A \code{character(1)}, either \code{"roc"} (ROC curves
#'   from the last CV fold) or \code{"bar"} (bar chart of mean CV AUC
#'   values). Default: \code{"roc"}.
#' @param colours A named \code{character} vector of colours for
#'   \code{"host_only"}, \code{"microbiome_only"}, and \code{"joint"}.
#'   Default: a colour-blind-friendly palette.
#' @param show_diagonal Logical. If \code{TRUE}, overlay the chance-level
#'   diagonal (AUC = 0.5 reference line). Default: \code{TRUE}.
#' @param linewidth A \code{numeric(1)} line width for ROC curves.
#'   Default: \code{0.9}.
#'
#' @return A \code{ggplot2} object.
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
#' result  <- MultiOmicsBridgeAnalysis(mae, outcome,
#'                                      n_components = 2,
#'                                      n_features_host = 20,
#'                                      n_features_mb = 10,
#'                                      n_biomarkers = 15, cv_folds = 3)
#' plotClassifierComparison(result, type = "bar")
#'
#' @seealso \code{\link{MultiOmicsBridgeAnalysis}},
#'   \code{\link{diagnosticClassifier}}, \code{\link{performance}}
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_abline geom_bar geom_errorbar
#' @importFrom ggplot2 geom_text scale_color_manual scale_fill_manual labs theme_bw theme
#' @importFrom ggplot2 element_text element_blank position_dodge coord_flip
#' @importFrom methods is
#' @importFrom stats setNames
#' @export
plotClassifierComparison <- function(result,
                                      type         = c("roc", "bar"),
                                      colours      = c(
                                          host_only       = "#4CAF50",
                                          microbiome_only = "#FF9800",
                                          joint           = "#2196F3"
                                      ),
                                      show_diagonal = TRUE,
                                      linewidth     = 0.9) {
    type <- match.arg(type)
    if (!is(result, "MOBResult"))
        stop("'result' must be a MOBResult object.")

    cr <- performance(result)
    if (length(cr) == 0L)
        stop("No classifier results found. Run MultiOmicsBridgeAnalysis() first.")

    model_names <- c("host_only", "microbiome_only", "joint")
    labels <- c(
        host_only       = "Host only",
        microbiome_only = "Microbiome only",
        joint           = "Joint (multi-omics)"
    )

    if (type == "bar") {
        # ── Bar chart of mean AUC ─────────────────────────────────────────────
        bar_df <- do.call(rbind, lapply(model_names, function(nm) {
            r <- cr[[nm]]
            if (is.null(r)) return(NULL)
            data.frame(
                model    = labels[nm],
                mean_auc = r$mean_auc %||% NA_real_,
                sd_auc   = r$sd_auc   %||% 0,
                stringsAsFactors = FALSE
            )
        }))
        bar_df <- bar_df[!is.na(bar_df$mean_auc), ]
        bar_df$model <- factor(bar_df$model, levels = rev(labels))

        fill_vals <- setNames(colours, labels[model_names])

        ggplot(bar_df,
               aes(x = .data[["model"]], y = .data[["mean_auc"]],
                   fill = .data[["model"]])) +
            geom_bar(stat = "identity", width = 0.55, colour = "white") +
            geom_errorbar(
                aes(ymin = .data[["mean_auc"]] - .data[["sd_auc"]],
                    ymax = .data[["mean_auc"]] + .data[["sd_auc"]]),
                width = 0.2, linewidth = 0.5, colour = "grey30"
            ) +
            geom_text(
                aes(label = sprintf("%.3f", .data[["mean_auc"]])),
                hjust = -0.15, size = 3.5
            ) +
            scale_fill_manual(values = fill_vals, guide = "none") +
            coord_flip() +
            labs(
                x     = NULL,
                y     = "Mean AUC-ROC (CV)",
                title = "Classifier Comparison: Multi-Omics Advantage",
                subtitle = sprintf("%d-fold cross-validation",
                                   cr$cv_folds %||% 5L)
            ) +
            theme_bw(base_size = 11) +
            theme(
                panel.grid.major.y = element_blank(),
                panel.grid.minor   = element_blank(),
                axis.text          = element_text(size = 10),
                plot.title         = element_text(size = 12, face = "bold"),
                plot.subtitle      = element_text(size = 10, colour = "grey40")
            )
    } else {
        # ── ROC curve overlay ─────────────────────────────────────────────────
        # Build ROC data frames from pROC roc objects (last CV fold)
        roc_dfs <- do.call(rbind, lapply(model_names, function(nm) {
            r <- cr[[nm]]
            if (is.null(r) || is.null(r$roc_data)) return(NULL)
            roc_obj <- r$roc_data
            data.frame(
                specificity = rev(roc_obj$specificities),
                sensitivity = rev(roc_obj$sensitivities),
                model       = unname(labels[nm]),
                auc         = as.numeric(pROC::auc(roc_obj)),
                stringsAsFactors = FALSE
            )
        }))

        if (is.null(roc_dfs) || nrow(roc_dfs) == 0L) {
            # Fallback: bar chart if no ROC data available
            message("ROC curve data not available. Plotting bar chart instead.")
            return(plotClassifierComparison(result, type = "bar",
                                             colours = colours))
        }

        roc_dfs$model <- factor(roc_dfs$model, levels = labels[model_names])
        # Add AUC to legend labels
        auc_labs <- vapply(model_names, function(nm) {
            r <- cr[[nm]]
            sprintf("%s (AUC=%.3f)", labels[nm], r$mean_auc %||% NA)
        }, character(1L))
        auc_labs_named <- setNames(auc_labs, labels[model_names])

        line_cols <- setNames(colours, labels[model_names])

        p <- ggplot(roc_dfs,
                    aes(x = 1 - .data[["specificity"]],
                        y = .data[["sensitivity"]],
                        colour = .data[["model"]])) +
            geom_line(linewidth = linewidth)

        if (show_diagonal) {
            p <- p + geom_abline(intercept = 0, slope = 1,
                                 linetype = "dashed",
                                 colour = "grey60", linewidth = 0.4)
        }

        p +
            scale_color_manual(values = line_cols, labels = auc_labs_named) +
            labs(
                x      = "1 - Specificity (FPR)",
                y      = "Sensitivity (TPR)",
                colour = "Classifier",
                title  = "ROC Curves: Multi-Omics Advantage",
                subtitle = sprintf("%d-fold cross-validation (last fold shown)",
                                   cr$cv_folds %||% 5L)
            ) +
            theme_bw(base_size = 11) +
            theme(
                panel.grid.minor = element_blank(),
                axis.text        = element_text(size = 10),
                legend.position  = "bottom",
                plot.title       = element_text(size = 12, face = "bold"),
                plot.subtitle    = element_text(size = 10, colour = "grey40")
            )
    }
}
