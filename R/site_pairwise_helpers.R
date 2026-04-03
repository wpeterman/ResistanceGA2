#' Build Pairwise Environmental-Difference Covariates
#'
#' Convert site-level environmental variables into pairwise covariates that are
#' aligned with the lower-triangle ordering used throughout \pkg{ResistanceGA2}.
#'
#' @param data Numeric vector, matrix, or data frame of site-level
#'   environmental variables.
#' @param vars Optional character or integer index selecting columns from
#'   \code{data}.
#' @param scale Logical. If \code{TRUE}, standardize each site-level variable
#'   before computing pairwise differences. Constant variables are converted to
#'   zeros.
#' @param transform Transformation applied to pairwise differences:
#'   \code{"absolute"} (default) or \code{"squared"}.
#' @param keep Optional logical/integer vector used to retain a subset of the
#'   final pairwise observations.
#' @param pop_n Optional integer vector giving the number of sampled
#'   individuals per population. When supplied, population-level environmental
#'   differences are expanded to the individual-level pair structure used by
#'   \code{\link[ResistanceGA2]{jl.prep}} with \code{pop2ind}.
#' @param suffix Optional suffix appended to output column names. Defaults to
#'   \code{"_diff"} for absolute differences and \code{"_sqdiff"} for squared
#'   differences.
#'
#' @return A data frame of pairwise environmental-difference covariates.
#'
#' @details
#' These helpers are intended for isolation-by-environment workflows, where
#' environmental mismatch between sample locations is used as one or more fixed
#' effects in an MLPE model. The output can be supplied directly to the
#' \code{covariates} argument of \code{\link[ResistanceGA2]{gdist.prep}} or
#' \code{\link[ResistanceGA2]{jl.prep}}.
#'
#' @examples
#' env <- data.frame(
#'   temp = c(10, 12, 15, 11),
#'   precip = c(100, 120, 80, 110)
#' )
#'
#' site_pairwise_diff(env)
#' site_pairwise_diff(env, vars = "temp", scale = TRUE)
#'
#' @export
site_pairwise_diff <- function(data,
                               vars = NULL,
                               scale = FALSE,
                               transform = c("absolute", "squared"),
                               keep = NULL,
                               pop_n = NULL,
                               suffix = NULL) {
  dat <- .site_pairwise_prepare_data(data, vars = vars)
  transform <- match.arg(transform)

  if (isTRUE(scale)) {
    dat <- .site_pairwise_scale_df(dat)
  }

  if (is.null(suffix)) {
    suffix <- if (identical(transform, "absolute")) "_diff" else "_sqdiff"
  }

  if (!is.character(suffix) || length(suffix) != 1L) {
    stop("'suffix' must be a single character string.")
  }

  out <- lapply(names(dat), function(nm) {
    mat <- outer(dat[[nm]], dat[[nm]], "-")

    if (identical(transform, "absolute")) {
      mat <- abs(mat)
    } else {
      mat <- mat ^ 2
    }

    .site_pairwise_finalize(mat, keep = keep, pop_n = pop_n)
  })

  out <- as.data.frame(out, stringsAsFactors = FALSE)
  names(out) <- paste0(names(dat), suffix)
  out
}


#' Build Pairwise Environmental-Distance Covariates
#'
#' Convert one or more site-level environmental variables into a single
#' pairwise environmental-distance covariate aligned with the lower-triangle
#' ordering used throughout \pkg{ResistanceGA2}.
#'
#' @param data Numeric vector, matrix, or data frame of site-level
#'   environmental variables.
#' @param vars Optional character or integer index selecting columns from
#'   \code{data}.
#' @param scale Logical. If \code{TRUE}, standardize each site-level variable
#'   before calculating environmental distance. Constant variables are
#'   converted to zeros.
#' @param method Distance metric passed to \code{\link[stats]{dist}}.
#' @param keep Optional logical/integer vector used to retain a subset of the
#'   final pairwise observations.
#' @param pop_n Optional integer vector giving the number of sampled
#'   individuals per population. When supplied, population-level environmental
#'   distances are expanded to the individual-level pair structure used by
#'   \code{\link[ResistanceGA2]{jl.prep}} with \code{pop2ind}.
#' @param name Name of the output covariate column. Default = \code{"env_dist"}.
#' @param ... Additional arguments passed to \code{\link[stats]{dist}}, such as
#'   \code{p} when \code{method = "minkowski"}.
#'
#' @return A one-column data frame containing pairwise environmental distance.
#'
#' @examples
#' env <- data.frame(
#'   temp = c(10, 12, 15, 11),
#'   precip = c(100, 120, 80, 110)
#' )
#'
#' site_pairwise_dist(env)
#' site_pairwise_dist(env, method = "manhattan", name = "climate_dist")
#'
#' @export
site_pairwise_dist <- function(data,
                               vars = NULL,
                               scale = TRUE,
                               method = "euclidean",
                               keep = NULL,
                               pop_n = NULL,
                               name = "env_dist",
                               ...) {
  dat <- .site_pairwise_prepare_data(data, vars = vars)

  if (isTRUE(scale)) {
    dat <- .site_pairwise_scale_df(dat)
  }

  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop("'name' must be a single, non-empty character string.")
  }

  mat <- as.matrix(stats::dist(dat, method = method, ...))
  out <- .site_pairwise_finalize(mat, keep = keep, pop_n = pop_n)

  stats::setNames(data.frame(out, stringsAsFactors = FALSE), name)
}


.site_pairwise_prepare_data <- function(data, vars = NULL) {
  if (is.atomic(data) && is.null(dim(data))) {
    data <- data.frame(x = as.numeric(data))
  } else if (is.matrix(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  } else if (!is.data.frame(data)) {
    stop("'data' must be a numeric vector, matrix, or data frame.")
  }

  if (!is.null(vars)) {
    data <- .site_pairwise_select_vars(data, vars = vars)
  }

  if (nrow(data) < 2L) {
    stop("'data' must contain at least two sites.")
  }

  if (ncol(data) < 1L) {
    stop("'data' must contain at least one environmental variable.")
  }

  if (!all(vapply(data, is.numeric, logical(1)))) {
    stop("All selected environmental variables must be numeric.")
  }

  if (anyNA(data)) {
    stop("Missing values are not allowed in site-level environmental data.")
  }

  data
}


.site_pairwise_select_vars <- function(data, vars) {
  if (is.character(vars)) {
    if (!all(vars %in% names(data))) {
      stop("All elements of 'vars' must match column names in 'data'.")
    }

    return(data[, vars, drop = FALSE])
  }

  if (is.numeric(vars)) {
    vars <- as.integer(vars)
    if (any(vars < 1L) || any(vars > ncol(data))) {
      stop("Numeric 'vars' indices are out of bounds for 'data'.")
    }

    return(data[, vars, drop = FALSE])
  }

  stop("'vars' must be NULL, a character vector of column names, or an integer vector of column indices.")
}


.site_pairwise_scale_df <- function(data) {
  out <- data

  for (nm in names(out)) {
    out[[nm]] <- .site_pairwise_scale_vec(out[[nm]])
  }

  out
}


.site_pairwise_scale_vec <- function(x) {
  s <- stats::sd(x)

  if (!is.finite(s) || s < sqrt(.Machine$double.eps)) {
    return(rep(0, length(x)))
  }

  as.numeric((x - mean(x)) / s)
}


.site_pairwise_finalize <- function(mat, keep = NULL, pop_n = NULL) {
  out <- if (is.null(pop_n)) lower(mat) else expand.mat(mat, pop_n)

  if (!is.null(keep)) {
    keep <- .mlpe_normalize_keep(keep, length(out), arg = "keep")
    out <- out[keep]
  }

  out
}
