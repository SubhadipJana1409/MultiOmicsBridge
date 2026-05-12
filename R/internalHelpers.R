# ── Internal helper functions for MultiOmicsBridge ────────────────────────────
# These functions are not exported. All exported API lives in the
# individual module files.

# ── Normalization helpers ─────────────────────────────────────────────────────

#' CLR transformation for microbiome count data
#'
#' Applies the centered log-ratio transformation to a taxa x samples
#' count matrix, adding a pseudocount to handle zero counts.
#'
#' @param x A \code{matrix} with taxa in rows and samples in columns.
#' @param pseudocount A positive \code{numeric} pseudocount added to all
#'   counts before log-transformation. Default: \code{0.5}.
#'
#' @return A \code{matrix} of CLR-transformed values (same dimensions as
#'   input).
#'
#' @keywords internal
#' @noRd
.clrTransform <- function(x, pseudocount = 0.5) {
  if (!is.matrix(x)) x <- as.matrix(x)
  x_pseudo <- x + pseudocount
  log_x <- log(x_pseudo)
  # Subtract column-wise geometric mean (mean of log values per sample)
  sweep(log_x, 2L, colMeans(log_x), FUN = "-")
}

#' Total sum scaling for microbiome count data
#'
#' Normalizes each sample to a fixed sum (proportions).
#'
#' @param x A \code{matrix} with taxa in rows and samples in columns.
#' @param scale_to Total sum per sample. Default: \code{1e6} (per-million).
#'
#' @return A \code{matrix} of TSS-normalized values.
#'
#' @keywords internal
#' @noRd
.tssTransform <- function(x, scale_to = 1e6) {
  if (!is.matrix(x)) x <- as.matrix(x)
  col_sums <- colSums(x)
  col_sums[col_sums == 0L] <- 1L
  sweep(x, 2L, col_sums / scale_to, FUN = "/")
}

#' TMM normalization + voom transformation for RNA-seq count data
#'
#' @param counts A \code{matrix} of raw counts, genes in rows and
#'   samples in columns.
#' @param design Optional design matrix for voom. Default: \code{NULL}.
#'
#' @return A list with \code{E} (voom log2-CPM matrix) and
#'   \code{weights} (observation-level precision weights).
#'
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom limma voom
#'
#' @keywords internal
#' @noRd
.voomNormalize <- function(counts, design = NULL) {
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  v <- limma::voom(dge, design = design, plot = FALSE)
  list(E = v$E, weights = v$weights)
}

# ── Cross-omics correlation helpers ──────────────────────────────────────────

#' Compute cross-omics Spearman correlation matrix
#'
#' Computes pairwise Spearman correlations between all host genes and all
#' microbial taxa across samples.
#'
#' @param host_mat A \code{matrix} (samples x genes).
#' @param mb_mat A \code{matrix} (samples x taxa).
#'
#' @return A \code{matrix} (genes x taxa) of Spearman correlation values.
#'
#' @keywords internal
#' @noRd
.crossOmicsCorrelation <- function(host_mat, mb_mat) {
  # Rank each column
  host_ranks <- apply(host_mat, 2L, rank, ties.method = "average")
  mb_ranks <- apply(mb_mat, 2L, rank, ties.method = "average")

  # Pearson correlation of ranks = Spearman correlation
  n <- nrow(host_mat)
  if (n < 3L) {
    return(matrix(NA_real_,
      nrow = ncol(host_mat),
      ncol = ncol(mb_mat)
    ))
  }

  host_scaled <- scale(host_ranks)
  mb_scaled <- scale(mb_ranks)

  cor_mat <- crossprod(host_scaled, mb_scaled) / (n - 1L)
  rownames(cor_mat) <- colnames(host_mat)
  colnames(cor_mat) <- colnames(mb_mat)
  cor_mat
}

# ── Classification helpers ────────────────────────────────────────────────────

#' Stratified k-fold cross-validation indices
#'
#' @param outcome A \code{factor} or character vector of outcomes.
#' @param k Number of folds. Default: \code{5}.
#' @param seed Random seed. Default: \code{42}.
#'
#' @return An \code{integer} vector of fold assignments (length = n).
#'
#' @keywords internal
#' @noRd
.stratifiedFolds <- function(outcome, k = 5L, seed = 42L) {
  outcome_f <- factor(outcome)
  fold_ids <- integer(length(outcome))
  for (lv in levels(outcome_f)) {
    idx <- which(outcome_f == lv)
    fids <- sample(rep(seq_len(k), length.out = length(idx)))
    fold_ids[idx] <- fids
  }
  fold_ids
}

#' Train a ranger Random Forest and return OOB AUC
#'
#' @param feature_mat A \code{matrix} (samples x features).
#' @param outcome A \code{factor}.
#' @param fold_ids An \code{integer} vector of fold assignments.
#' @param seed Integer random seed. Default: \code{42}.
#'
#' @return A \code{list} with \code{mean_auc}, \code{sd_auc},
#'   \code{fold_auc}, and \code{roc_data} (for the last fold).
#'
#' @importFrom ranger ranger
#' @importFrom pROC roc auc
#' @importFrom stats predict sd
#'
#' @keywords internal
#' @noRd
.cvRandomForest <- function(feature_mat, outcome, fold_ids, seed = 42L) {
  k <- max(fold_ids)
  auc_vals <- numeric(k)
  roc_last <- NULL
  outcome_f <- factor(outcome)
  pos_class <- levels(outcome_f)[2L]

  for (f in seq_len(k)) {
    train_idx <- fold_ids != f
    test_idx <- fold_ids == f

    if (sum(train_idx) < 4L || sum(test_idx) < 2L) {
      auc_vals[f] <- NA_real_
      next
    }

    train_df <- as.data.frame(feature_mat[train_idx, , drop = FALSE])
    test_df <- as.data.frame(feature_mat[test_idx, , drop = FALSE])
    train_df[[".__outcome__."]] <- outcome_f[train_idx]

    rf_fit <- ranger::ranger(
      formula = .__outcome__. ~ .,
      data = train_df,
      num.trees = 500L,
      probability = TRUE,
      seed = seed + f,
      verbose = FALSE
    )

    pred <- predict(rf_fit, data = test_df)$predictions
    pred_p <- pred[, pos_class]
    true_y <- outcome_f[test_idx]

    roc_obj <- pROC::roc(
      response = true_y,
      predictor = pred_p,
      levels = levels(outcome_f),
      direction = "<",
      quiet = TRUE
    )
    auc_vals[f] <- as.numeric(pROC::auc(roc_obj))
    if (f == k) roc_last <- roc_obj
  }

  auc_vals_clean <- auc_vals[!is.na(auc_vals)]
  list(
    mean_auc = mean(auc_vals_clean),
    sd_auc   = if (length(auc_vals_clean) > 1L) sd(auc_vals_clean) else 0,
    fold_auc = auc_vals,
    roc_data = roc_last
  )
}

# ── DIABLO helpers ────────────────────────────────────────────────────────────

#' Extract sample scores from DIABLO result
#'
#' @param diablo_res Result object from \code{mixOmics::block.splsda}.
#' @param which_block Name of the block to extract scores from.
#'   Use \code{"AVE_inner"} for the consensus. Default: \code{"host"}.
#'
#' @return A \code{matrix} (samples x components) of scores.
#'
#' @keywords internal
#' @noRd
.extractScores <- function(diablo_res, which_block = "host") {
  variates <- diablo_res[["variates"]]
  if (which_block %in% names(variates)) {
    variates[[which_block]]
  } else {
    variates[[1L]]
  }
}

#' Extract feature loadings from DIABLO result
#'
#' @param diablo_res Result object from \code{mixOmics::block.splsda}.
#' @param block_name Name of the block ("host" or "microbiome").
#'
#' @return A \code{matrix} (features x components) of loadings.
#'
#' @keywords internal
#' @noRd
.extractLoadings <- function(diablo_res, block_name) {
  loadings <- diablo_res[["loadings"]]
  if (block_name %in% names(loadings)) {
    loadings[[block_name]]
  } else {
    matrix(nrow = 0L, ncol = 0L)
  }
}

#' Build ranked biomarker DataFrame from loadings
#'
#' @param host_loadings A \code{matrix} (genes x components).
#' @param mb_loadings A \code{matrix} (taxa x components).
#' @param n_biomarkers Integer, top N per layer.
#'
#' @return A \code{DataFrame} with biomarker rankings.
#'
#' @importFrom S4Vectors DataFrame
#'
#' @keywords internal
#' @noRd
.rankBiomarkers <- function(host_loadings, mb_loadings, n_biomarkers = 50L) {
  .rank_layer <- function(loadings, layer_name) {
    if (nrow(loadings) == 0L || ncol(loadings) == 0L) {
      return(data.frame())
    }

    # Combined loading score: sqrt sum of squares across components
    scores <- sqrt(rowSums(loadings^2))
    top_n <- min(n_biomarkers, length(scores))
    ord <- order(scores, decreasing = TRUE)[seq_len(top_n)]

    # Find which component has the highest absolute loading per gene
    best_comp <- apply(abs(loadings), 1L, which.max)

    data.frame(
      feature = rownames(loadings)[ord],
      omics_layer = layer_name,
      loading_score = scores[ord],
      rank = seq_len(top_n),
      component = best_comp[ord],
      stringsAsFactors = FALSE
    )
  }

  host_df <- .rank_layer(host_loadings, "host")
  mb_df <- .rank_layer(mb_loadings, "microbiome")
  combined <- rbind(host_df, mb_df)

  if (nrow(combined) == 0L) {
    return(S4Vectors::DataFrame())
  }

  S4Vectors::DataFrame(combined)
}
