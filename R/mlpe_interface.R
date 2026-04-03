#' Build Long-Form Data for MLPE Models
#'
#' Convert pairwise responses and covariates into a long-form data frame that is
#' ready for \code{\link[ResistanceGA2]{mlpe}}. Inputs can be supplied either as
#' symmetric matrices/\code{dist} objects, or as long vectors paired with an
#' explicit two-column pair index.
#'
#' @param response Pairwise response as a lower-triangle vector, symmetric
#'   matrix, \code{dist} object, or a long vector aligned with \code{pairs}.
#' @param ... Named pairwise covariates supplied in the same format as
#'   \code{response}.
#' @param pairs Optional data frame or matrix containing at least two columns
#'   that identify the endpoints for each observation. When supplied,
#'   \code{response} and all covariates are treated as long vectors aligned to
#'   these rows. Additional columns in \code{pairs} are preserved.
#' @param labels Optional vector of endpoint labels used when converting a
#'   symmetric matrix, \code{dist} object, or lower-triangle vector.
#' @param pair_names Names to use for the generated endpoint columns when
#'   \code{pairs} is not supplied.
#' @param keep Optional logical/integer vector used to retain a subset of pair
#'   observations.
#'
#' @return A data frame with one row per retained pairwise observation.
#'
#' @importFrom methods as
#'
#' @examples
#' y <- matrix(rnorm(25), 5)
#' x <- matrix(rnorm(25), 5)
#'
#' dat <- mlpe_data(response = y, resistance = x)
#'
#' head(dat)
#'
#' id <- data.frame(from = c(1, 1, 2), to = c(2, 3, 3))
#' mlpe_data(response = c(0.2, 0.4, 0.6), x = c(1, 2, 3), pairs = id)
#'
#' @export
mlpe_data <- function(response,
                      ...,
                      pairs = NULL,
                      labels = NULL,
                      pair_names = c("from", "to"),
                      keep = NULL) {
  covariates <- list(...)

  if (!is.null(pairs)) {
    pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)

    if (ncol(pairs) < 2) {
      stop("'pairs' must contain at least two columns.")
    }

    response_vec <- .mlpe_long_vector(response, expected_n = nrow(pairs), arg = "response")
    data <- pairs

    for (nm in names(covariates)) {
      data[[nm]] <- .mlpe_long_vector(covariates[[nm]],
                                      expected_n = nrow(pairs),
                                      arg = nm)
    }
  } else {
    pair_index <- .mlpe_pair_index(response,
                                   labels = labels,
                                   pair_names = pair_names,
                                   arg = "response")
    response_vec <- pair_index$values
    data <- pair_index$pairs

    for (nm in names(covariates)) {
      data[[nm]] <- .mlpe_pair_values(covariates[[nm]],
                                      labels = pair_index$labels,
                                      expected_n = length(response_vec),
                                      arg = nm)
    }
  }

  data$response <- response_vec

  if (!is.null(keep)) {
    keep <- .mlpe_normalize_keep(keep, nrow(data), arg = "keep")
    data <- data[keep, , drop = FALSE]
  }

  rownames(data) <- NULL
  data
}


#' Fit MLPE Models With an lme4-Style Formula
#'
#' Fit maximum-likelihood population effects (MLPE) models using standard
#' \code{lme4}-style formulas while defining dyadic random-intercept terms from
#' explicit endpoint columns.
#'
#' @param formula Mixed-effects model formula.
#' @param data Data frame containing the response, covariates, and endpoint
#'   columns.
#' @param pairs Definition of one or more MLPE dyad terms. Supply either a
#'   character vector of length two, which is interpreted as
#'   \code{list(pair = c("from", "to"))}, or a named list mapping random-effect
#'   term names in \code{formula} to the two endpoint columns used to construct
#'   each dyadic random intercept.
#' @param REML Logical. If TRUE, fit by restricted maximum likelihood for linear
#'   mixed models.
#' @param keep Optional logical/integer vector used to retain a subset of rows
#'   from \code{data} before fitting.
#' @param ... Additional arguments passed to \code{\link[lme4]{glFormula}} when
#'   fitting generalized models, such as \code{family}.
#'
#' @return A fitted \code{merMod} object.
#'
#' @details
#' Dyadic random effects are specified in \code{formula} using standard random
#' intercept syntax, for example \code{(1 | pair)}. The corresponding endpoint
#' columns are declared through \code{pairs}, e.g.
#' \code{pairs = c("from", "to")} or
#' \code{pairs = list(pair = c("from", "to"))}.
#'
#' Additional ordinary \code{lme4} random effects may be included in the same
#' model. The dyadic covariance contribution is rebuilt after any row subsetting
#' so omitted pairwise observations are handled correctly.
#'
#' @examples
#' y <- matrix(rnorm(25), 5)
#' x <- matrix(rnorm(25), 5)
#'
#' dat <- mlpe_data(response = y, resistance = x)
#'
#' fit <- mlpe(
#'   response ~ resistance + (1 | pair),
#'   data = dat,
#'   pairs = c("from", "to")
#' )
#'
#' @export
mlpe <- function(formula,
                 data,
                 pairs,
                 REML = FALSE,
                 keep = NULL,
                 ...) {
  .mlpe_fit_mermod(formula = formula,
                   data = data,
                   REML = REML,
                   keep = keep,
                   pairs = pairs,
                   ...)
}


.mlpe_attach_workflow_pairs <- function(data, ID) {
  if (is.null(data) || is.null(ID)) {
    return(data)
  }

  if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  if (nrow(data) != nrow(ID)) {
    stop("Workflow MLPE metadata requires 'data' and 'ID' to have the same number of rows.")
  }

  pair_terms <- list()

  if (all(c("pop1.ind", "pop2.ind", "pop1.pop", "pop2.pop") %in% names(ID))) {
    data$.mlpe_ind1 <- ID$pop1.ind
    data$.mlpe_ind2 <- ID$pop2.ind
    data$.mlpe_grp1 <- ID$pop1.pop
    data$.mlpe_grp2 <- ID$pop2.pop

    pair_terms$pop <- c(".mlpe_ind1", ".mlpe_ind2")
    pair_terms$grp <- c(".mlpe_grp1", ".mlpe_grp2")

    if ("corr_" %in% names(ID)) {
      data$.mlpe_corr <- ID$corr_
      pair_terms$cor.grp <- list(
        cols = c(".mlpe_ind1", ".mlpe_ind2"),
        active = ".mlpe_corr"
      )
    }
  } else if (all(c("pop1", "pop2") %in% names(ID))) {
    data$.mlpe_pop1 <- ID$pop1
    data$.mlpe_pop2 <- ID$pop2

    pair_terms$pop <- c(".mlpe_pop1", ".mlpe_pop2")

    if ("corr_" %in% names(ID)) {
      data$.mlpe_corr <- ID$corr_
      pair_terms$cor.grp <- list(
        cols = c(".mlpe_pop1", ".mlpe_pop2"),
        active = ".mlpe_corr"
      )
    }
  }

  if (length(pair_terms) > 0L) {
    attr(data, "mlpe_pairs") <- pair_terms
  }

  data
}


.mlpe_formula_from_pair_terms <- function(response, predictor, pair_terms) {
  rhs <- c(
    predictor,
    sprintf("(1 | %s)", names(pair_terms))
  )

  as.formula(paste(response, "~", paste(rhs, collapse = " + ")))
}


.mlpe_formula_from_data <- function(data,
                                    response,
                                    predictor,
                                    fallback = NULL) {
  pair_terms <- attr(data, "mlpe_pairs", exact = TRUE)

  if (is.list(pair_terms) && length(pair_terms) > 0L) {
    return(.mlpe_formula_from_pair_terms(response, predictor, pair_terms))
  }

  if (is.null(fallback)) {
    stop("No MLPE pair metadata is attached to 'data', and no fallback formula was supplied.")
  }

  fallback
}


.mlpe_fit_mermod <- function(formula,
                             data,
                             REML = FALSE,
                             keep = NULL,
                             pairs = NULL,
                             ZZ = NULL,
                             ...) {
  if (!inherits(formula, "formula")) {
    formula <- as.formula(formula)
  }

  if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  n_before_keep <- nrow(data)

  if (!is.null(keep)) {
    keep <- .mlpe_normalize_keep(keep, n_before_keep, arg = "keep")
    data <- data[keep, , drop = FALSE]
  }

  data <- .mlpe_drop_unused_factors(data)

  if (!is.null(pairs) && !is.null(ZZ)) {
    stop("Use either 'pairs' or 'ZZ', not both.")
  }

  pair_terms <- .mlpe_normalize_pairs(pairs, data)
  if (length(pair_terms) > 0L) {
    data <- .mlpe_add_placeholders(data, pair_terms)
  }

  args <- list(...)
  generalized <- any("family" %in% names(args))
  mod <- .mlpe_make_mod(formula = formula,
                        data = data,
                        REML = REML,
                        generalized = generalized,
                        args = args)

  if (length(pair_terms) > 0L) {
    mod$reTrms$Zt <- .mlpe_replace_pair_blocks(mod$reTrms,
                                               data = data,
                                               pair_terms = pair_terms)
  } else if (!is.null(ZZ)) {
    ZZ <- .mlpe_subset_ZZ_columns(ZZ = ZZ,
                                  keep = keep,
                                  expected_ncol = nrow(data))

    if (nrow(mod$reTrms$Zt) != nrow(ZZ)) {
      data <- .mlpe_expand_legacy_placeholder(data = data,
                                              reTrms = mod$reTrms,
                                              target_nrow = nrow(ZZ))
      mod <- .mlpe_make_mod(formula = formula,
                            data = data,
                            REML = REML,
                            generalized = generalized,
                            args = args)
    }

    mod$reTrms$Zt <- .mlpe_align_ZZ(ZZ = ZZ,
                                    expected_ncol = nrow(data),
                                    expected_nrow = nrow(mod$reTrms$Zt))
  }

  if (generalized) {
    dfun <- do.call(lme4::mkGlmerDevfun, mod)
    opt <- lme4::optimizeGlmer(dfun)
  } else {
    dfun <- do.call(lme4::mkLmerDevfun, mod)
    opt <- lme4::optimizeLmer(dfun)
  }

  lme4::mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr)
}


.mlpe_make_mod <- function(formula, data, REML, generalized, args) {
  if (generalized) {
    if (isTRUE(REML)) {
      cat("REML will be ignored when fitting a generalized MLPE model.")
    }

    do.call(lme4::glFormula, c(list(formula = formula, data = data), args))
  } else {
    lme4::lFormula(formula,
                   data = data,
                   REML = REML)
  }
}


.mlpe_replace_pair_blocks <- function(reTrms, data, pair_terms) {
  term_names <- names(reTrms$cnms)
  blocks <- vector("list", length(term_names))

  for (i in seq_along(term_names)) {
    nm <- term_names[[i]]
    rows <- seq.int(reTrms$Gp[[i]] + 1L, reTrms$Gp[[i + 1L]])

    if (nm %in% names(pair_terms)) {
      if (!identical(reTrms$cnms[[nm]], "(Intercept)")) {
        stop("MLPE dyad terms must be random intercepts only. Problem term: '",
             nm, "'.")
      }

      levels_nm <- levels(reTrms$flist[[nm]])
      blocks[[i]] <- .mlpe_pair_block(data = data,
                                      pair_term = pair_terms[[nm]],
                                      levels = levels_nm,
                                      term = nm)
    } else {
      blocks[[i]] <- reTrms$Zt[rows, , drop = FALSE]
    }
  }

  out <- do.call(rbind, blocks)
  Matrix::drop0(out)
}


.mlpe_pair_block <- function(data, pair_term, levels, term) {
  pair_cols <- unname(pair_term$cols)
  active <- .mlpe_pair_active(data, pair_term, term)
  x1 <- data[[pair_cols[[1]]]]
  x2 <- data[[pair_cols[[2]]]]

  if (any(is.na(x1)) || any(is.na(x2))) {
    stop("Missing endpoint labels are not allowed for MLPE term '", term, "'.")
  }

  x1_chr <- as.character(x1)
  x2_chr <- as.character(x2)

  if (any(x1_chr == x2_chr)) {
    stop("Self-comparisons are not allowed for MLPE term '", term, "'.")
  }

  f1 <- factor(x1_chr, levels = levels)
  f2 <- factor(x2_chr, levels = levels)

  if (any(is.na(f1)) || any(is.na(f2))) {
    stop("Endpoint labels for MLPE term '", term,
         "' are not consistent with the model term levels.")
  }

  Z1 <- Matrix::fac2sparse(f1, "d", drop = FALSE)
  Z2 <- Matrix::fac2sparse(f2, "d", drop = FALSE)

  if (any(!active)) {
    Z1[, !active] <- 0
    Z2[, !active] <- 0
  }

  Matrix::drop0(Z1 + Z2)
}


.mlpe_expand_legacy_placeholder <- function(data, reTrms, target_nrow) {
  term_names <- names(reTrms$cnms)

  if (length(term_names) != 1L || !identical(reTrms$cnms[[1]], "(Intercept)")) {
    stop("Custom 'ZZ' does not align with the random-effects structure implied by the formula.")
  }

  if (nrow(data) <= target_nrow) {
    stop("Custom 'ZZ' requires ", target_nrow,
         " sampled grouping levels, but only ", nrow(data),
         " observations remain after filtering. The legacy lme4-backed path requires more observations than grouping levels.")
  }

  data[[term_names[[1]]]] <- factor(rep(seq_len(target_nrow), length.out = nrow(data)))
  data
}


.mlpe_subset_ZZ_columns <- function(ZZ, keep, expected_ncol) {
  if (!inherits(ZZ, "Matrix")) {
    ZZ <- methods::as(ZZ, "dgCMatrix")
  }

  if (!is.null(keep) && ncol(ZZ) != expected_ncol) {
    if (length(keep) == ncol(ZZ)) {
      ZZ <- ZZ[, keep, drop = FALSE]
    }
  }

  if (ncol(ZZ) != expected_ncol) {
    stop("Custom 'ZZ' has ", ncol(ZZ), " columns, but the fitted model expects ",
         expected_ncol, ".")
  }

  ZZ
}


.mlpe_align_ZZ <- function(ZZ, keep = NULL, expected_ncol, expected_nrow) {
  ZZ <- .mlpe_subset_ZZ_columns(ZZ = ZZ,
                                keep = keep,
                                expected_ncol = expected_ncol)
  if (nrow(ZZ) != expected_nrow) {
    keep_rows <- Matrix::rowSums(ZZ != 0) > 0
    ZZ <- ZZ[keep_rows, , drop = FALSE]
  }

  if (nrow(ZZ) != expected_nrow) {
    stop("Custom 'ZZ' has ", nrow(ZZ), " rows after alignment, but the fitted model expects ",
         expected_nrow, ".")
  }

  Matrix::drop0(ZZ)
}


.mlpe_normalize_pairs <- function(pairs, data) {
  if (is.null(pairs)) {
    return(list())
  }

  if (is.character(pairs) && length(pairs) == 2L) {
    pairs <- list(pair = pairs)
  }

  if (!is.list(pairs) || is.null(names(pairs)) || any(names(pairs) == "")) {
    stop("'pairs' must be NULL, a character vector of length two, or a named list of endpoint-column pairs.")
  }

  out <- lapply(names(pairs), function(nm) {
    pair_term <- pairs[[nm]]

    if (is.character(pair_term) && length(pair_term) == 2L) {
      pair_term <- list(cols = pair_term)
    }

    if (!is.list(pair_term) || is.null(pair_term$cols)) {
      stop("Each element of 'pairs' must be a character vector of length two or a list with a 'cols' element.")
    }

    cols <- pair_term$cols
    if (!is.character(cols) || length(cols) != 2L) {
      stop("Each element of 'pairs' must define exactly two endpoint columns.")
    }

    if (!all(cols %in% names(data))) {
      stop("Endpoint columns for MLPE term '", nm, "' were not found in 'data'.")
    }

    active <- pair_term$active
    if (!is.null(active)) {
      if (!is.character(active) || length(active) != 1L) {
        stop("The 'active' field for MLPE term '", nm, "' must name a single column in 'data'.")
      }

      if (!active %in% names(data)) {
        stop("Active-row column '", active, "' for MLPE term '", nm, "' was not found in 'data'.")
      }
    }

    list(cols = cols, active = active)
  })

  names(out) <- names(pairs)
  out
}


.mlpe_add_placeholders <- function(data, pair_terms) {
  for (nm in names(pair_terms)) {
    pair_term <- pair_terms[[nm]]
    cols <- pair_term$cols
    active <- .mlpe_pair_active(data, pair_term, nm)
    levels_nm <- .mlpe_pair_levels(data[[cols[[1]]]],
                                   data[[cols[[2]]]],
                                   active = active)

    if (length(levels_nm) < 2L) {
      stop("MLPE term '", nm, "' must involve at least two distinct endpoints.")
    }

    if (nrow(data) <= length(levels_nm)) {
      stop(
        "MLPE term '", nm, "' has ", length(levels_nm),
        " endpoint levels but only ", nrow(data),
        " observations. The lme4-backed interface requires more observations than endpoint levels after filtering."
      )
    }

    existing_levels <- NULL
    if (nm %in% names(data)) {
      existing_levels <- unique(as.character(data[[nm]]))
      existing_levels <- existing_levels[!is.na(existing_levels)]
    }

    if (!identical(sort(existing_levels), sort(levels_nm))) {
      data[[nm]] <- factor(rep(levels_nm, length.out = nrow(data)),
                           levels = levels_nm)
    }
  }

  data
}


.mlpe_pair_levels <- function(x1, x2, active) {
  used_levels <- unique(c(as.character(x1[active]), as.character(x2[active])))

  if (is.factor(x1) || is.factor(x2)) {
    level_order <- unique(c(
      if (is.factor(x1)) levels(x1) else as.character(x1[active]),
      if (is.factor(x2)) levels(x2) else as.character(x2[active])
    ))

    return(level_order[level_order %in% used_levels])
  }

  used_levels
}


.mlpe_pair_active <- function(data, pair_term, term) {
  if (is.null(pair_term$active)) {
    return(rep(TRUE, nrow(data)))
  }

  active <- as.logical(data[[pair_term$active]])
  if (length(active) != nrow(data) || any(is.na(active))) {
    stop("Active-row column for MLPE term '", term, "' must be a non-missing logical or 0/1 vector with one value per observation.")
  }

  active
}


.mlpe_drop_unused_factors <- function(data) {
  for (nm in names(data)) {
    if (is.factor(data[[nm]])) {
      data[[nm]] <- droplevels(data[[nm]])
    }
  }

  data
}


.mlpe_normalize_keep <- function(keep, n, arg = "keep") {
  if (length(keep) != n) {
    stop("'", arg, "' must have length ", n, ".")
  }

  keep <- as.logical(keep)

  if (any(is.na(keep))) {
    stop("'", arg, "' cannot contain NA values.")
  }

  keep
}


.mlpe_pair_index <- function(x, labels = NULL, pair_names = c("from", "to"), arg = "response") {
  if (inherits(x, "dist")) {
    n <- attr(x, "Size")
    x_labels <- attr(x, "Labels")
    values <- as.vector(x)
  } else if (is.matrix(x) || is.data.frame(x)) {
    x <- as.matrix(x)
    if (nrow(x) != ncol(x)) {
      stop("'", arg, "' must be a square matrix when supplied as a matrix.")
    }

    n <- nrow(x)
    x_labels <- rownames(x)
    values <- lower(x)
  } else if (is.vector(x)) {
    values <- as.vector(x)
    n <- .mlpe_infer_nodes(length(values), arg = arg)
    x_labels <- NULL
  } else {
    stop("Unsupported input type for '", arg, "'.")
  }

  if (is.null(labels)) {
    labels <- x_labels
  }

  if (is.null(labels)) {
    labels <- seq_len(n)
  }

  if (length(labels) != n) {
    stop("'labels' must have length ", n, ".")
  }

  idx <- .mlpe_lower_index(n)
  pairs <- data.frame(
    labels[idx$from],
    labels[idx$to],
    stringsAsFactors = FALSE
  )
  names(pairs)[1:2] <- pair_names

  list(values = values, pairs = pairs, labels = labels)
}


.mlpe_pair_values <- function(x, labels, expected_n, arg = "x") {
  out <- .mlpe_pair_index(x, labels = labels, arg = arg)

  if (length(out$values) != expected_n) {
    stop("'", arg, "' does not align with the response pair structure.")
  }

  out$values
}


.mlpe_long_vector <- function(x, expected_n, arg = "x") {
  if (is.matrix(x) || is.data.frame(x) || inherits(x, "dist")) {
    stop("'", arg, "' must be a long vector when 'pairs' is supplied.")
  }

  x <- as.vector(x)
  if (length(x) != expected_n) {
    stop("'", arg, "' must have length ", expected_n, ".")
  }

  x
}


.mlpe_infer_nodes <- function(n_pairs, arg = "response") {
  n <- (1 + sqrt(1 + 8 * n_pairs)) / 2

  if (!isTRUE(all.equal(n, round(n)))) {
    stop("Length of '", arg, "' does not correspond to the lower triangle of a symmetric matrix.")
  }

  as.integer(round(n))
}


.mlpe_lower_index <- function(n) {
  mat <- matrix(0, nrow = n, ncol = n)
  list(
    from = col(mat)[lower.tri(mat)],
    to = row(mat)[lower.tri(mat)]
  )
}
