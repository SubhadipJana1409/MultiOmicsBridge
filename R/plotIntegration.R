#' @title Joint Biplot of Multi-Omics Integration Results
#'
#' @description
#' Produces a scatter plot of integrated sample scores in the space
#' defined by two DIABLO latent components, with samples coloured by
#' outcome group and optional feature loading vectors overlaid as arrows.
#' This biplot makes it immediately clear how well the multi-omics
#' integration separates the outcome groups and which features drive
#' that separation.
#'
#' @param result A \code{\link{MOBResult}} object from
#'   \code{\link{MultiOmicsBridgeAnalysis}}.
#' @param comp A length-2 \code{integer} vector specifying which two
#'   components to plot. Default: \code{c(1, 2)}.
#' @param outcome A \code{character} or \code{factor} vector of outcome
#'   labels (one per sample), used for point colouring. If \code{NULL},
#'   points are plotted in a single colour.
#' @param show_loadings Logical. If \code{TRUE}, overlay the top feature
#'   loading vectors as arrows. Default: \code{TRUE}.
#' @param n_loading_arrows An \code{integer(1)} number of top loading
#'   arrows to display per omics layer. Default: \code{5}.
#' @param point_size A \code{numeric(1)} point size. Default: \code{2.5}.
#' @param point_alpha A \code{numeric(1)} point transparency. Default:
#'   \code{0.8}.
#' @param colours A named \code{character} vector of colours for the
#'   outcome groups. If \code{NULL}, a default palette is used.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' library(S4Vectors)
#' scores <- matrix(rnorm(40), nrow = 20, ncol = 2,
#'     dimnames = list(paste0("S", 1:20), c("Comp1","Comp2")))
#' fl <- list(
#'     host = matrix(rnorm(20), nrow = 10, ncol = 2,
#'         dimnames = list(paste0("G", 1:10), c("Comp1","Comp2"))),
#'     microbiome = matrix(rnorm(10), nrow = 5, ncol = 2,
#'         dimnames = list(paste0("T", 1:5), c("Comp1","Comp2")))
#' )
#' bm <- DataFrame(feature = c("G1","T1"), omics_layer = c("host","microbiome"),
#'                 loading_score = c(0.8,0.6), rank = c(1L,2L), component = c(1L,1L))
#' obj <- MOBResult(scores, fl, bm, list(),
#'                  params = list(outcome_levels = c("ctrl","treat")))
#' outcome <- rep(c("ctrl","treat"), each = 10)
#' plotIntegration(obj, outcome = outcome)
#'
#' @seealso \code{\link{MultiOmicsBridgeAnalysis}}, \code{\link{MOBResult}},
#'   \code{\link{plotClassifierComparison}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_text
#' @importFrom ggplot2 scale_color_manual scale_color_brewer labs theme_bw theme
#' @importFrom ggplot2 element_text element_blank geom_hline geom_vline
#' @importFrom methods is
#' @importFrom utils head
#' @importFrom stats setNames
#' @export
plotIntegration <- function(result,
                             comp             = c(1L, 2L),
                             outcome          = NULL,
                             show_loadings    = TRUE,
                             n_loading_arrows = 5L,
                             point_size       = 2.5,
                             point_alpha      = 0.8,
                             colours          = NULL) {
    if (!is(result, "MOBResult"))
        stop("'result' must be a MOBResult object.")

    scores <- integrationScores(result)
    if (is.null(scores) || nrow(scores) == 0L)
        stop("No integration scores found in 'result'.")
    if (ncol(scores) < max(comp))
        stop("'comp' requests component ", max(comp),
             " but only ", ncol(scores), " available.")

    c1 <- comp[1L]; c2 <- comp[2L]
    cn1 <- colnames(scores)[c1]
    cn2 <- colnames(scores)[c2]

    df <- data.frame(
        x       = scores[, c1],
        y       = scores[, c2],
        sample  = rownames(scores),
        outcome = if (!is.null(outcome)) as.character(outcome) else "all",
        stringsAsFactors = FALSE
    )

    # Default colours
    if (is.null(colours)) {
        n_groups <- length(unique(df$outcome))
        colours  <- setNames(
            c("#2196F3","#F44336","#4CAF50","#FF9800")[seq_len(n_groups)],
            unique(df$outcome)
        )
    }

    p <- ggplot(df, aes(x = .data[["x"]], y = .data[["y"]],
                        colour = .data[["outcome"]])) +
        geom_hline(yintercept = 0, linetype = "dashed",
                   colour = "grey70", linewidth = 0.3) +
        geom_vline(xintercept = 0, linetype = "dashed",
                   colour = "grey70", linewidth = 0.3) +
        geom_point(size = point_size, alpha = point_alpha) +
        scale_color_manual(values = colours) +
        labs(
            x      = paste0("Component ", c1),
            y      = paste0("Component ", c2),
            colour = "Group",
            title  = "Multi-Omics Joint Integration (DIABLO)",
            subtitle = if (!is.null(result@params$integration_method))
                paste0("Method: ", result@params$integration_method) else NULL
        ) +
        theme_bw(base_size = 11) +
        theme(
            panel.grid.minor = element_blank(),
            axis.text        = element_text(size = 10),
            plot.title       = element_text(size = 12, face = "bold")
        )

    # Overlay loading arrows
    if (show_loadings && n_loading_arrows > 0L) {
        fl <- featureLoadings(result)
        arrow_df <- do.call(rbind, lapply(names(fl), function(layer) {
            lmat <- fl[[layer]]
            if (is.null(lmat) || nrow(lmat) == 0L || ncol(lmat) < max(comp))
                return(NULL)
            scores_l2 <- sqrt(rowSums(lmat^2))
            top_idx   <- order(scores_l2, decreasing = TRUE)[
                seq_len(min(n_loading_arrows, nrow(lmat)))]
            scale_f   <- max(abs(scores)) / max(abs(lmat[top_idx, ]))
            clean_labels <- function(x) {
                x <- gsub("^[a-z]__", "", x)
                ifelse(nchar(x) > 20, paste0(substr(x, 1, 18), "..."), x)
            }
            data.frame(
                xend  = lmat[top_idx, c1] * scale_f * 0.8,
                yend  = lmat[top_idx, c2] * scale_f * 0.8,
                label = clean_labels(rownames(lmat)[top_idx]),
                layer = layer,
                stringsAsFactors = FALSE
            )
        }))

        if (!is.null(arrow_df) && nrow(arrow_df) > 0L) {
            arrow_df$x <- 0; arrow_df$y <- 0
            
            # Colour arrows by layer (Host = Blue, Microbiome = Orange)
            arrow_df$arr_col <- ifelse(arrow_df$layer == "host", "#2196F3", "#FF9800")
            
            p <- p +
                geom_segment(
                    data = arrow_df,
                    aes(x = .data[["x"]], y = .data[["y"]],
                        xend = .data[["xend"]], yend = .data[["yend"]]),
                    arrow = grid::arrow(length = grid::unit(0.15, "cm")),
                    colour = arrow_df$arr_col, linewidth = 0.5,
                    inherit.aes = FALSE
                ) +
                ggrepel::geom_text_repel(
                    data = arrow_df,
                    aes(x = .data[["xend"]], y = .data[["yend"]],
                        label = .data[["label"]]),
                    size = 3.0, colour = arrow_df$arr_col, fontface = "bold",
                    segment.colour = NA, box.padding = 0.5, force = 2,
                    inherit.aes = FALSE
                )
        }
    }
    p
}
