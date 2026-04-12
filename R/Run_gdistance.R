#' Run gdistance cost/commute distance calculation
#'
#' Calculates pairwise cost or commute distances across a resistance surface
#' using the \code{gdistance} package.
#'
#' @param gdist.inputs Object created by \code{\link[ResistanceGA2]{gdist.prep}}.
#' @param r A single-layer \code{SpatRaster} object (from the \pkg{terra}
#'   package).
#' @param scl Logical. Scale the correction values (default \code{TRUE}). Set
#'   to \code{FALSE} to obtain absolute distance values. See
#'   \code{\link[gdistance]{geoCorrection}} for details.
#' @param return.error.value Logical. If \code{TRUE}, return \code{-99999} on
#'   warnings/errors instead of throwing an error. This is intended for
#'   internal GA fitness evaluation, where failed candidate surfaces should be
#'   penalized but optimization should continue. Default is \code{FALSE} for
#'   user-facing calls.
#' @param commute.approx Optional override for the commute-distance
#'   approximation stored in \code{gdist.inputs}. \code{NULL} uses the setting
#'   from \code{\link{gdist.prep}}. \code{"none"} runs the exact gdistance
#'   commute solve; \code{"aggregate"} coarsens \code{r} first for faster
#'   iterative optimization.
#' @param approx.factor Optional integer aggregation factor override used when
#'   \code{commute.approx = "aggregate"}.
#' @param approx.scale Optional logical override. If \code{TRUE}, aggregate
#'   approximation distances are multiplied by \code{approx.factor^2}.
#'
#' @return A numeric distance vector (or matrix) of pairwise cost/commute
#'   distances. If \code{return.error.value = TRUE}, returns \code{-99999} on
#'   error or warning.
#'
#' @details
#' \pkg{gdistance} still operates on \pkg{raster} objects. This function keeps
#' the public API terra-based by coercing the supplied \code{SpatRaster}
#' internally before building the transition object. The default transition
#' function uses vectorized sparse builders for projected 4- and 8-neighbor
#' transition matrices and projected geographic corrections. Selected-pair runs
#' reuse a single distance calculation over the required endpoints.
#' Commute-distance runs solve the reduced Laplacian for all requested endpoints
#' in one batch.
#'
#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @examples
#' pts <- terra::vect(samples[, 2:3], type = "points")
#' gdist.inputs <- gdist.prep(
#'   n.Pops = nrow(samples),
#'   samples = pts,
#'   method = "costDistance"
#' )
#'
#' r_surface <- terra::subset(terra::unwrap(resistance_surfaces), "continuous")
#' cd <- Run_gdistance(gdist.inputs, r_surface)
#' length(cd)
Run_gdistance <- function(gdist.inputs,
                          r,
                          scl = TRUE,
                          return.error.value = FALSE,
                          commute.approx = NULL,
                          approx.factor = NULL,
                          approx.scale = NULL) {
  out <- tryCatch(
    {
      approx <- .gdistance_commute_approx_settings(
        gdist.inputs,
        commute.approx = commute.approx,
        approx.factor = approx.factor,
        approx.scale = approx.scale
      )

      if (identical(gdist.inputs$method, "commuteDistance") &&
          identical(approx$method, "aggregate") &&
          approx$factor > 1L) {
        r <- .gdistance_aggregate_raster(r, approx$factor)
      }

      r_gd <- .gdistance_raster(r, arg = "r")

      fast_projected <- .gdistance_use_fast_mean_transition(r_gd, gdist.inputs)
      tr <- .gdistance_transition(
        r_gd,
        gdist.inputs,
        correction = fast_projected,
        scl = scl
      )

      samples <- gdist.inputs$samples  # coordinate matrix
      if (is.null(gdist.inputs$keep)) {
        pair_map <- list(pairs = NULL)
        run_samples <- samples
      } else {
        pair_map <- .gdistance_pair_map(gdist.inputs)
        run_samples <- samples[pair_map$endpoints, , drop = FALSE]
        if (!length(pair_map$endpoints)) {
          return(numeric(0))
        }
      }

      if (gdist.inputs$method == 'costDistance') {
        trC <- if (fast_projected) tr else .gdistance_geoCorrection(tr, "c", scl = scl)
        ret <- gdistance::costDistance(trC, run_samples)
        rm(trC)
      } else {
        trR <- if (fast_projected) tr else .gdistance_geoCorrection(tr, "r", scl = scl)
        ret <- .gdistance_commuteDistance_fast(trR, run_samples) / 1000
        if (identical(approx$method, "aggregate") && isTRUE(approx$scale)) {
          ret <- ret * (approx$factor^2)
        }
        rm(trR)
      }

      if (!is.null(pair_map$pairs)) {
        ret <- .gdistance_extract_pairs(ret, pair_map)
      }

      rm(tr, r_gd)
      gc()
      return(ret)
    },
    warning = function(w) {
      return(-99999)
    },
    error = function(e) {
      return(-99999)
    }
  )
  if (identical(out, -99999) && !isTRUE(return.error.value)) {
    stop("Run_gdistance failed for the supplied `gdist.inputs` and resistance surface.")
  }
  out
}

.gdistance_commute_approx_settings <- function(gdist.inputs,
                                               commute.approx = NULL,
                                               approx.factor = NULL,
                                               approx.scale = NULL) {
  method <- commute.approx
  if (is.null(method)) {
    method <- gdist.inputs$commute.approx
  }
  if (is.null(method)) {
    method <- "none"
  }
  method <- match.arg(method, c("none", "aggregate"))

  factor <- approx.factor
  if (is.null(factor)) {
    factor <- gdist.inputs$approx.factor
  }
  if (is.null(factor)) {
    factor <- 4L
  }
  factor <- as.integer(factor)
  if (length(factor) != 1L || is.na(factor) || factor < 1L) {
    stop("'approx.factor' must be a single positive integer.")
  }

  scale <- approx.scale
  if (is.null(scale)) {
    scale <- gdist.inputs$approx.scale
  }
  if (is.null(scale)) {
    scale <- TRUE
  }

  list(method = method, factor = factor, scale = isTRUE(scale))
}

.gdistance_aggregate_raster <- function(r, factor) {
  r <- .validate_spatraster(r, arg = "r", nlyr = 1L)
  terra::aggregate(r, fact = factor, fun = mean, na.rm = TRUE)
}

.gdistance_transition <- function(r_gd, gdist.inputs, correction = FALSE, scl = TRUE) {
  if (.gdistance_use_fast_mean_transition(r_gd, gdist.inputs)) {
    return(.gdistance_transition_mean(
      r_gd,
      gdist.inputs$directions,
      correction = correction,
      scl = scl
    ))
  }

  gdistance::transition(
    x = r_gd,
    transitionFunction = gdist.inputs$transitionFunction,
    directions = gdist.inputs$directions
  )
}

.gdistance_transition_cache <- new.env(parent = emptyenv())
.gdistance_transition_cache$keys <- character()

.gdistance_use_fast_mean_transition <- function(r_gd, gdist.inputs) {
  isTRUE(gdist.inputs$directions %in% c(4, 8)) &&
    isFALSE(raster::isLonLat(r_gd)) &&
    is.function(gdist.inputs$transitionFunction) &&
    identical(body(gdist.inputs$transitionFunction), body(function(x) 1 / mean(x))) &&
    identical(names(formals(gdist.inputs$transitionFunction)), "x")
}

.gdistance_transition_mean <- function(r_gd, directions, correction = FALSE, scl = TRUE) {
  template <- .gdistance_transition_template(r_gd, directions)
  values <- raster::getValues(r_gd)
  keep <- is.finite(values[template$from]) & is.finite(values[template$to])

  edge_from <- template$from[keep]
  edge_to <- template$to[keep]
  edge_values <- 2 / (values[edge_from] + values[edge_to])
  keep_values <- edge_values != 0

  edge_from <- edge_from[keep_values]
  edge_to <- edge_to[keep_values]
  edge_values <- edge_values[keep_values]

  if (isTRUE(correction)) {
    correction_values <- if (isTRUE(scl)) {
      template$geo_factor_scaled[keep][keep_values]
    } else {
      template$geo_factor_unscaled[keep][keep_values]
    }
    edge_values <- edge_values * correction_values
  }

  if (!all(edge_values >= 0)) {
    warning("transition function gives negative values")
  }

  transition_matrix <- Matrix::sparseMatrix(
    i = edge_from,
    j = edge_to,
    x = edge_values,
    dims = c(template$ncells, template$ncells),
    symmetric = TRUE
  )

  methods::new(
    "TransitionLayer",
    nrows = as.integer(template$nrows),
    ncols = as.integer(template$ncols),
    extent = raster::extent(r_gd),
    crs = raster::projection(r_gd, asText = FALSE),
    transitionMatrix = transition_matrix,
    transitionCells = seq_len(template$ncells),
    matrixValues = "conductance"
  )
}

.gdistance_transition_template <- function(r_gd, directions) {
  key <- .gdistance_transition_template_key(r_gd, directions)
  if (exists(key, envir = .gdistance_transition_cache, inherits = FALSE)) {
    return(get(key, envir = .gdistance_transition_cache, inherits = FALSE))
  }

  nrows <- raster::nrow(r_gd)
  ncols <- raster::ncol(r_gd)
  ncells <- raster::ncell(r_gd)
  cells <- seq_len(ncells)
  rows <- ((cells - 1L) %/% ncols) + 1L
  cols <- ((cells - 1L) %% ncols) + 1L

  if (directions == 4) {
    offsets <- c(-ncols, -1L)
    valid <- list(rows > 1L, cols > 1L)
  } else {
    offsets <- c(-ncols - 1L, -ncols, -ncols + 1L, -1L)
    valid <- list(
      rows > 1L & cols > 1L,
      rows > 1L,
      rows > 1L & cols < ncols,
      cols > 1L
    )
  }

  edge_from <- vector("list", length(offsets))
  edge_to <- vector("list", length(offsets))

  for (i in seq_along(offsets)) {
    to_cells <- cells[valid[[i]]]
    from_cells <- to_cells + offsets[[i]]
    edge_from[[i]] <- from_cells
    edge_to[[i]] <- to_cells
  }

  edge_from <- unlist(edge_from, use.names = FALSE)
  edge_to <- unlist(edge_to, use.names = FALSE)
  edge_order <- order(edge_to, edge_from)
  edge_from <- edge_from[edge_order]
  edge_to <- edge_to[edge_order]

  row1 <- ((edge_from - 1L) %/% ncols) + 1L
  row2 <- ((edge_to - 1L) %/% ncols) + 1L
  col1 <- ((edge_from - 1L) %% ncols) + 1L
  col2 <- ((edge_to - 1L) %% ncols) + 1L
  edge_distance <- sqrt(
    ((col1 - col2) * raster::xres(r_gd))^2 +
      ((row1 - row2) * raster::yres(r_gd))^2
  )

  midpoint <- c(
    mean(c(raster::xmin(r_gd), raster::xmax(r_gd))),
    mean(c(raster::ymin(r_gd), raster::ymax(r_gd)))
  )
  scale_value <- raster::pointDistance(
    midpoint,
    midpoint + c(raster::xres(r_gd), 0),
    longlat = FALSE
  )

  template <- list(
    nrows = nrows,
    ncols = ncols,
    ncells = ncells,
    from = edge_from,
    to = edge_to,
    geo_factor_scaled = scale_value / edge_distance,
    geo_factor_unscaled = 1 / edge_distance
  )

  .gdistance_transition_cache_set(key, template)
  template
}

.gdistance_transition_template_key <- function(r_gd, directions) {
  paste(
    "mean",
    directions,
    raster::nrow(r_gd),
    raster::ncol(r_gd),
    format(raster::xmin(r_gd), digits = 17),
    format(raster::xmax(r_gd), digits = 17),
    format(raster::ymin(r_gd), digits = 17),
    format(raster::ymax(r_gd), digits = 17),
    raster::projection(r_gd),
    sep = "|"
  )
}

.gdistance_transition_cache_set <- function(key, template) {
  assign(key, template, envir = .gdistance_transition_cache)
  keys <- .gdistance_transition_cache$keys
  if (!key %in% keys) {
    keys <- c(keys, key)
  }
  while (length(keys) > 3L) {
    rm(list = keys[[1]], envir = .gdistance_transition_cache)
    keys <- keys[-1L]
  }
  .gdistance_transition_cache$keys <- keys
}

.gdistance_geoCorrection <- function(x, type, scl = TRUE) {
  if (.gdistance_use_fast_geoCorrection(x, type)) {
    return(.gdistance_geoCorrection_projected(x, scl = scl))
  }

  gdistance::geoCorrection(x, type, scl = scl)
}

.gdistance_use_fast_geoCorrection <- function(x, type) {
  type %in% c("c", "r") &&
    isFALSE(raster::isLonLat(x)) &&
    identical(gdistance::matrixValues(x), "conductance") &&
    methods::is(gdistance::transitionMatrix(x), "CsparseMatrix")
}

.gdistance_geoCorrection_projected <- function(x, scl = TRUE) {
  transition_matrix <- gdistance::transitionMatrix(x)

  scale_value <- 1
  if (isTRUE(scl)) {
    midpoint <- c(
      mean(c(raster::xmin(x), raster::xmax(x))),
      mean(c(raster::ymin(x), raster::ymax(x)))
    )
    scale_value <- raster::pointDistance(
      midpoint,
      midpoint + c(raster::xres(x), 0),
      longlat = FALSE
    )
  }

  col_counts <- diff(transition_matrix@p)
  col_cells <- rep.int(seq_len(transition_matrix@Dim[2]), col_counts)
  row_cells <- transition_matrix@i + 1L

  row1 <- ((row_cells - 1L) %/% raster::ncol(x)) + 1L
  row2 <- ((col_cells - 1L) %/% raster::ncol(x)) + 1L
  col1 <- ((row_cells - 1L) %% raster::ncol(x)) + 1L
  col2 <- ((col_cells - 1L) %% raster::ncol(x)) + 1L

  edge_distance <- sqrt(
    ((col1 - col2) * raster::xres(x))^2 +
      ((row1 - row2) * raster::yres(x))^2
  )

  transition_matrix@x <- transition_matrix@x * (scale_value / edge_distance)
  gdistance::transitionMatrix(x) <- transition_matrix
  x
}

.gdistance_pair_map <- function(gdist.inputs) {
  n <- nrow(gdist.inputs$samples)

  if (is.null(gdist.inputs$keep)) {
    return(list(endpoints = seq_len(n), pairs = NULL))
  }

  keep <- which(gdist.inputs$keep == 1)
  if (!length(keep)) {
    return(list(endpoints = integer(0), pairs = matrix(integer(0), ncol = 2)))
  }

  id <- as.matrix(gdist.inputs$ID[, 1:2, drop = FALSE])
  id <- apply(id, 2, function(x) as.integer(as.character(x)))
  id <- id[keep, , drop = FALSE]

  endpoints <- sort(unique(as.integer(id)))
  pairs <- cbind(
    match(id[, 1], endpoints),
    match(id[, 2], endpoints)
  )

  list(endpoints = endpoints, pairs = pairs)
}

.gdistance_extract_pairs <- function(distance, pair_map) {
  if (!length(pair_map$pairs)) {
    return(numeric(0))
  }

  distance <- as.matrix(distance)
  as.numeric(distance[pair_map$pairs])
}

.gdistance_commuteDistance_fast <- function(x, coords) {
  x_orig <- x
  transition_matrix <- gdistance::transitionMatrix(x)
  if (isFALSE(methods::is(transition_matrix, "dsCMatrix"))) {
    stop("symmetric transition matrix required",
         "(dsCMatrix) in TransitionLayer object x")
  }

  coords <- as.matrix(coords)
  rd <- matrix(NA_real_, nrow = nrow(coords), ncol = nrow(coords))
  rownames(rd) <- rownames(coords)
  colnames(rd) <- rownames(coords)

  all_from_cells <- raster::cellFromXY(x, coords)
  if (!all(!is.na(all_from_cells))) {
    return(gdistance::commuteDistance(x_orig, coords))
  }

  x <- .gdistance_transition_solidify(x)
  from_cells <- all_from_cells[all_from_cells %in% gdistance::transitionCells(x)]
  if (length(from_cells) < length(all_from_cells)) {
    return(gdistance::commuteDistance(x_orig, coords))
  }

  from_cells <- unique(all_from_cells)
  laplacian <- .gdistance_laplacian(x)
  n <- max(laplacian@Dim)
  laplacian_reduced <- laplacian[-n, -n]

  index <- match(from_cells, gdistance::transitionCells(x))
  if (anyNA(index) || any(index >= n)) {
    return(gdistance::commuteDistance(x_orig, coords))
  }

  Lplus <- .gdistance_lplus_at_cells(laplacian_reduced, index, n)

  diag_Lplus <- diag(Lplus)
  rdSS <- outer(diag_Lplus, diag_Lplus, "+") - (2 * Lplus)
  rdSS <- rdSS * sum(gdistance::transitionMatrix(x))

  index1 <- which(all_from_cells %in% from_cells)
  index2 <- match(all_from_cells[all_from_cells %in% from_cells], from_cells)
  rd[index1, index1] <- rdSS[index2, index2]

  rd <- stats::as.dist(rd)
  attr(rd, "method") <- "commute"
  rd
}

.gdistance_lplus_at_cells <- function(laplacian_reduced, index, n) {
  nr <- nrow(laplacian_reduced)
  C <- 1e-300 * n
  k <- length(index)
  block_size <- min(k, 64L)
  Lplus <- matrix(NA_real_, nrow = k, ncol = k)

  chol <- tryCatch(
    Matrix::Cholesky(
      laplacian_reduced,
      LDL = FALSE,
      super = nr >= 250000L
    ),
    error = function(e) NULL
  )

  for (start in seq.int(1L, k, by = block_size)) {
    cols <- seq.int(start, min(start + block_size - 1L, k))
    rhs <- matrix(-C / n, nrow = nr, ncol = length(cols))
    rhs[cbind(index[cols], seq_along(cols))] <- C - (C / n)

    xi <- if (is.null(chol)) {
      Matrix::solve(laplacian_reduced, rhs)
    } else {
      Matrix::solve(chol, rhs)
    }

    endpoint_rows <- as.matrix(xi[index, , drop = FALSE])
    Lplus[, cols] <- sweep(
      endpoint_rows,
      2,
      Matrix::colSums(xi) / n,
      FUN = "-"
    ) / C
  }

  Lplus
}

.gdistance_transition_solidify <- function(x) {
  transition_matrix <- gdistance::transitionMatrix(x, inflate = FALSE)
  selection <- which(Matrix::rowMeans(transition_matrix) > 1e-300)
  x@transitionCells <- gdistance::transitionCells(x)[selection]
  x@transitionMatrix <- transition_matrix[selection, selection]
  x
}

.gdistance_laplacian <- function(x) {
  transition_matrix <- gdistance::transitionMatrix(x, inflate = FALSE)
  laplacian <- Matrix::Diagonal(x = Matrix::colSums(transition_matrix)) -
    transition_matrix
  methods::as(laplacian, "symmetricMatrix")
}
