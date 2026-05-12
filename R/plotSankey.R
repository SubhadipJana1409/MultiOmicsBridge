#' @title Feature Flow Diagram: Omics Layer to Biomarkers to Outcome
#'
#' @description
#' Produces a Sankey-style flow diagram showing the pipeline from data
#' source (host or microbiome) through selected top biomarkers to
#' predicted outcome classes. The width of each connection is proportional
#' to the feature's loading score, making it easy to see which features
#' contribute most to separating the outcome groups. The diagram uses
#' base \code{ggplot2} geometry (no additional Sankey packages required).
#'
#' @param result A \code{\link{MOBResult}} object from
#'   \code{\link{MultiOmicsBridgeAnalysis}}.
#' @param n_features An \code{integer(1)} number of top features to
#'   display per omics layer. Default: \code{10}.
#' @param colours A named \code{character} vector with colours for
#'   \code{"host"}, \code{"microbiome"}, and each outcome level.
#'   If \code{NULL}, defaults are used.
#' @param node_width A \code{numeric(1)} width of the node rectangles.
#'   Default: \code{0.15}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' library(S4Vectors)
#' set.seed(42)
#' scores <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'     dimnames = list(paste0("S", 1:10), c("Comp1","Comp2")))
#' fl <- list(
#'     host = matrix(abs(rnorm(20)), nrow = 10, ncol = 2,
#'         dimnames = list(paste0("Gene", 1:10), c("Comp1","Comp2"))),
#'     microbiome = matrix(abs(rnorm(10)), nrow = 5, ncol = 2,
#'         dimnames = list(paste0("Taxon", 1:5), c("Comp1","Comp2")))
#' )
#' bm <- DataFrame(
#'     feature = c(paste0("Gene", 1:5), paste0("Taxon", 1:3)),
#'     omics_layer = c(rep("host", 5), rep("microbiome", 3)),
#'     loading_score = c(0.9, 0.7, 0.5, 0.4, 0.3, 0.8, 0.6, 0.2),
#'     rank = 1:8, component = rep(1L, 8)
#' )
#' obj <- MOBResult(scores, fl, bm, list(),
#'     params = list(outcome_levels = c("ctrl","treat")))
#' plotSankey(obj, n_features = 5)
#'
#' @seealso \code{\link{MultiOmicsBridgeAnalysis}},
#'   \code{\link{plotIntegration}}, \code{\link{plotBiomarkerNetwork}}
#'
#' @importFrom ggplot2 ggplot aes geom_rect geom_segment geom_text
#' @importFrom ggplot2 scale_fill_manual scale_colour_manual labs theme_void theme
#' @importFrom ggplot2 element_text
#' @importFrom methods is
#' @importFrom utils head
#' @importFrom stats setNames
#' @export
plotSankey <- function(result,
                        n_features = 10L,
                        colours    = NULL,
                        node_width = 0.15) {
    if (!is(result, "MOBResult"))
        stop("'result' must be a MOBResult object.")

    bm <- as.data.frame(biomarkers(result))
    if (nrow(bm) == 0L)
        stop("No biomarkers found in 'result'.")

    outcome_lvls <- result@params$outcome_levels %||% c("Group1", "Group2")

    # ── Select top features per layer ─────────────────────────────────────────
    host_bm <- bm[bm$omics_layer == "host", ]
    mb_bm   <- bm[bm$omics_layer == "microbiome", ]
    host_bm <- host_bm[order(host_bm$loading_score, decreasing = TRUE), ]
    mb_bm   <- mb_bm[order(mb_bm$loading_score, decreasing = TRUE), ]

    host_top <- utils::head(host_bm, n_features)
    mb_top   <- utils::head(mb_bm,   n_features)
    all_bm   <- rbind(host_top, mb_top)

    if (nrow(all_bm) == 0L) stop("No features to plot.")

    # -- Layout (3 columns: source -> biomarker -> outcome) --------------------
    # X positions: source = 0, biomarker = 1, outcome = 2
    x_src <- 0.05; x_bm <- 0.5; x_out <- 0.95

    # Source nodes (host, microbiome)
    src_nodes <- data.frame(
        label  = c("Host\nRNA-seq", "Microbiome"),
        x      = x_src,
        y      = c(0.75, 0.25),
        fill   = c("host", "microbiome"),
        stringsAsFactors = FALSE
    )

    # Biomarker nodes (evenly spaced vertically)
    n_bm  <- nrow(all_bm)
    bm_y  <- seq(0.95, 0.05, length.out = n_bm)
    clean_labels <- function(x) {
        x <- gsub("^[a-z]__", "", x)
        ifelse(nchar(x) > 22, paste0(substr(x, 1, 20), "..."), x)
    }

    bm_nodes <- data.frame(
        label   = clean_labels(all_bm$feature),
        x       = x_bm,
        y       = bm_y,
        fill    = all_bm$omics_layer,
        score   = all_bm$loading_score,
        stringsAsFactors = FALSE
    )

    # Outcome nodes
    out_nodes <- data.frame(
        label = outcome_lvls,
        x     = x_out,
        y     = c(0.7, 0.3),
        fill  = outcome_lvls,
        stringsAsFactors = FALSE
    )

    # ── Default colours ───────────────────────────────────────────────────────
    if (is.null(colours)) {
        colours <- c(
            host          = "#2196F3",
            microbiome    = "#FF9800",
            setNames(c("#F44336", "#4CAF50"), outcome_lvls)
        )
    }

    # -- Edges: source -> biomarker --------------------------------------------
    src_to_bm <- do.call(rbind, lapply(seq_len(n_bm), function(i) {
        src_y <- if (all_bm$omics_layer[i] == "host") 0.75 else 0.25
        data.frame(
            x    = x_src + node_width / 2,
            xend = x_bm  - node_width / 2,
            y    = src_y,
            yend = bm_y[i],
            lwd  = pmax(0.3, all_bm$loading_score[i] * 2),
            fill = all_bm$omics_layer[i],
            stringsAsFactors = FALSE
        )
    }))

    # -- Edges: biomarker -> outcome (each bm connects to both outcomes) -------
    bm_to_out <- do.call(rbind, lapply(seq_len(n_bm), function(i) {
        out_y <- c(0.7, 0.3)
        do.call(rbind, lapply(seq_along(outcome_lvls), function(j) {
            data.frame(
                x    = x_bm  + node_width / 2,
                xend = x_out - node_width / 2,
                y    = bm_y[i],
                yend = out_y[j],
                lwd  = pmax(0.2, all_bm$loading_score[i] * 1.2),
                fill = outcome_lvls[j],
                stringsAsFactors = FALSE
            )
        }))
    }))

    # ── All nodes combined ────────────────────────────────────────────────────
    all_nodes <- rbind(
        data.frame(label = src_nodes$label, x = src_nodes$x,
                   y = src_nodes$y, fill = src_nodes$fill,
                   stringsAsFactors = FALSE),
        data.frame(label = bm_nodes$label, x = bm_nodes$x,
                   y = bm_nodes$y, fill = bm_nodes$fill,
                   stringsAsFactors = FALSE),
        data.frame(label = out_nodes$label, x = out_nodes$x,
                   y = out_nodes$y, fill = out_nodes$fill,
                   stringsAsFactors = FALSE)
    )

    hw <- node_width / 2
    all_nodes$xmin <- all_nodes$x - hw
    all_nodes$xmax <- all_nodes$x + hw
    all_nodes$ymin <- all_nodes$y - 0.04
    all_nodes$ymax <- all_nodes$y + 0.04

    all_edges <- rbind(src_to_bm, bm_to_out)

    ggplot() +
        # Edges
        geom_segment(
            data = all_edges,
            aes(x = .data[["x"]], y = .data[["y"]],
                xend = .data[["xend"]], yend = .data[["yend"]],
                linewidth = .data[["lwd"]],
                colour = .data[["fill"]]),
            alpha = 0.25,
            lineend = "round"
        ) +
        # Nodes
        geom_rect(
            data = all_nodes,
            aes(xmin = .data[["xmin"]], xmax = .data[["xmax"]],
                ymin = .data[["ymin"]], ymax = .data[["ymax"]],
                fill = .data[["fill"]]),
            colour = "white", linewidth = 0.4
        ) +
        # Labels
        geom_text(
            data = all_nodes,
            aes(x = .data[["x"]], y = .data[["y"]],
                label = .data[["label"]]),
            size = 2.5, fontface = "bold", colour = "white"
        ) +
        scale_fill_manual(values = colours, guide = "none") +
        scale_colour_manual(values = colours, guide = "none") +
        labs(
            title    = "Feature Flow: Omics Layer -> Biomarkers -> Outcome",
            subtitle = sprintf("Top %d features per layer", n_features)
        ) +
        theme_void() +
        theme(
            plot.title    = element_text(size = 12, face = "bold",
                                         hjust = 0.5),
            plot.subtitle = element_text(size = 10, colour = "grey40",
                                         hjust = 0.5),
            plot.margin   = ggplot2::margin(10, 10, 10, 10)
        )
}
