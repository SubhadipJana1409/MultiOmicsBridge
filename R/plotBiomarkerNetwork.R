#' @title Cross-Omics Biomarker Correlation Heatmap
#'
#' @description
#' Produces a clustered heatmap of Spearman correlations between the top
#' selected host genes and microbial taxa. Rows represent host genes,
#' columns represent microbial taxa, and cell colour encodes the Spearman
#' correlation value. Strong positive or negative correlations between a
#' host gene and a microbial taxon suggest a potential functional
#' host-microbe interaction.
#'
#' @param result A \code{\link{MOBResult}} object from
#'   \code{\link{MultiOmicsBridgeAnalysis}}.
#' @param mae A \code{MultiAssayExperiment} from \code{\link{matchSamples}},
#'   used to compute correlations from assay data.
#' @param n_host An \code{integer(1)} number of top host genes to display.
#'   Default: \code{20}.
#' @param n_mb An \code{integer(1)} number of top microbial taxa to
#'   display. Default: \code{15}.
#' @param host_assay A \code{character(1)} assay name in host SE.
#'   Default: \code{"voom"}.
#' @param mb_assay A \code{character(1)} assay name in microbiome SE.
#'   Default: \code{"CLR"}.
#' @param cor_thresh A \code{numeric(1)} minimum absolute correlation to
#'   display (cells below this threshold are shown in grey). Default:
#'   \code{0}.
#' @param low_colour A \code{character(1)} colour for strong negative
#'   correlations. Default: \code{"#F44336"}.
#' @param mid_colour A \code{character(1)} colour for zero correlation.
#'   Default: \code{"white"}.
#' @param high_colour A \code{character(1)} colour for strong positive
#'   correlations. Default: \code{"#2196F3"}.
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
#' plotBiomarkerNetwork(result, mae, n_host = 10, n_mb = 8)
#'
#' @seealso \code{\link{MultiOmicsBridgeAnalysis}}, \code{\link{plotIntegration}},
#'   \code{\link{plotTCRatio}}
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 labs
#'   theme_bw theme element_text element_blank
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom methods is
#' @importFrom stats hclust dist
#' @export
plotBiomarkerNetwork <- function(result,
                                  mae,
                                  n_host      = 20L,
                                  n_mb        = 15L,
                                  host_assay  = "voom",
                                  mb_assay    = "CLR",
                                  cor_thresh  = 0,
                                  low_colour  = "#F44336",
                                  mid_colour  = "white",
                                  high_colour = "#2196F3") {
    if (!is(result, "MOBResult"))
        stop("'result' must be a MOBResult object.")
    if (!is(mae, "MultiAssayExperiment"))
        stop("'mae' must be a MultiAssayExperiment.")

    # ── Select features ───────────────────────────────────────────────────────
    bm <- as.data.frame(biomarkers(result))

    host_bm <- bm[bm$omics_layer == "host", ]
    host_bm <- host_bm[order(host_bm$loading_score, decreasing = TRUE), ]
    mb_bm   <- bm[bm$omics_layer == "microbiome", ]
    mb_bm   <- mb_bm[order(mb_bm$loading_score, decreasing = TRUE), ]

    host_sel <- utils::head(host_bm$feature, n_host)
    mb_sel   <- utils::head(mb_bm$feature, n_mb)

    if (length(host_sel) == 0L || length(mb_sel) == 0L)
        stop("No biomarkers found. Run MultiOmicsBridgeAnalysis() first.")

    # ── Extract assay matrices ────────────────────────────────────────────────
    host_se <- MultiAssayExperiment::experiments(mae)[["host"]]
    mb_se   <- MultiAssayExperiment::experiments(mae)[["microbiome"]]

    h_assay <- if (host_assay %in% SummarizedExperiment::assayNames(host_se))
        host_assay else SummarizedExperiment::assayNames(host_se)[1L]
    m_assay <- if (mb_assay %in% SummarizedExperiment::assayNames(mb_se))
        mb_assay else SummarizedExperiment::assayNames(mb_se)[1L]

    host_mat <- t(SummarizedExperiment::assay(host_se, h_assay))
    mb_mat   <- t(SummarizedExperiment::assay(mb_se,   m_assay))

    host_use <- intersect(host_sel, colnames(host_mat))
    mb_use   <- intersect(mb_sel,   colnames(mb_mat))

    if (length(host_use) < 2L || length(mb_use) < 2L)
        stop("Too few features after matching with assay matrices.")

    host_sub <- host_mat[, host_use, drop = FALSE]
    mb_sub   <- mb_mat[,   mb_use,   drop = FALSE]

    # ── Compute correlation matrix ────────────────────────────────────────────
    cor_mat <- .crossOmicsCorrelation(host_sub, mb_sub)

    # Threshold
    if (cor_thresh > 0)
        cor_mat[abs(cor_mat) < cor_thresh] <- NA_real_

    # Hierarchical clustering for row/col ordering
    if (nrow(cor_mat) > 2L && !all(is.na(cor_mat))) {
        hc_row  <- tryCatch(
            hclust(dist(cor_mat, method = "euclidean")),
            error = function(e) NULL)
        if (!is.null(hc_row))
            cor_mat <- cor_mat[hc_row$order, , drop = FALSE]
    }
    if (ncol(cor_mat) > 2L && !all(is.na(cor_mat))) {
        hc_col  <- tryCatch(
            hclust(dist(t(cor_mat), method = "euclidean")),
            error = function(e) NULL)
        if (!is.null(hc_col))
            cor_mat <- cor_mat[, hc_col$order, drop = FALSE]
    }

    # ── Melt to long format ───────────────────────────────────────────────────
    plot_df <- do.call(rbind, lapply(seq_len(nrow(cor_mat)), function(i) {
        data.frame(
            gene       = rownames(cor_mat)[i],
            taxon      = colnames(cor_mat),
            correlation = cor_mat[i, ],
            stringsAsFactors = FALSE
        )
    }))
    plot_df$gene  <- factor(plot_df$gene,  levels = rownames(cor_mat))
    plot_df$taxon <- factor(plot_df$taxon, levels = colnames(cor_mat))

    ggplot(plot_df, aes(x = .data[["taxon"]], y = .data[["gene"]],
                        fill = .data[["correlation"]])) +
        geom_tile(colour = "white", linewidth = 0.2) +
        scale_fill_gradient2(
            low      = low_colour,
            mid      = mid_colour,
            high     = high_colour,
            midpoint = 0,
            limits   = c(-1, 1),
            na.value = "grey85",
            name     = "Spearman r"
        ) +
        labs(
            x     = "Microbial taxon",
            y     = "Host gene",
            title = "Cross-Omics Biomarker Correlation Network",
            subtitle = sprintf("Top %d genes × %d taxa",
                               length(host_use), length(mb_use))
        ) +
        theme_bw(base_size = 10) +
        theme(
            axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y      = element_text(size = 8),
            panel.grid       = element_blank(),
            panel.border     = element_blank(),
            plot.title       = element_text(size = 12, face = "bold"),
            plot.subtitle    = element_text(size = 10, colour = "grey40")
        )
}
