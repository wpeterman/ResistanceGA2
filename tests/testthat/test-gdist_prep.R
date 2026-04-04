# Helper: make a coordinate matrix and response vector for n populations
make_coords <- function(n, seed = 1) {
  set.seed(seed)
  matrix(runif(n * 2, 0, 100), ncol = 2)
}

make_points <- function(n, seed = 1) {
  terra::vect(make_coords(n, seed = seed), type = "points")
}

make_gd <- function(n, seed = 2) {
  set.seed(seed)
  runif(choose(n, 2))
}

# ----- Basic structure --------------------------------------------------------

test_that("gdist.prep returns list with all expected names", {
  coords <- make_points(10)
  gd     <- make_gd(10)
  out    <- suppressWarnings(
    gdist.prep(n.Pops = 10, response = gd, samples = coords)
  )
  expected <- c("response", "samples", "covariates", "formula",
                "transitionFunction", "directions", "ID", "ZZ",
                "keep", "n.Pops", "longlat", "method", "df")
  expect_named(out, expected)
})

test_that("gdist.prep stores samples as matrix", {
  coords <- make_points(8)
  out    <- gdist.prep(n.Pops = 8, samples = coords)
  expect_true(is.matrix(out$samples))
  expect_equal(ncol(out$samples), 2L)
})

test_that("gdist.prep ZZ has dimensions n × choose(n,2)", {
  n      <- 7
  coords <- make_points(n)
  gd     <- make_gd(n)
  out    <- suppressWarnings(
    gdist.prep(n.Pops = n, response = gd, samples = coords)
  )
  expect_equal(nrow(out$ZZ), n)
  expect_equal(ncol(out$ZZ), choose(n, 2))
})

test_that("gdist.prep default formula is gd ~ cd + (1 | pop)", {
  coords <- make_points(10)
  gd     <- make_gd(10)
  out    <- suppressWarnings(
    gdist.prep(n.Pops = 10, response = gd, samples = coords)
  )
  expect_equal(deparse(out$formula), "gd ~ cd + (1 | pop)")
})

test_that("gdist.prep df has gd and pop columns when response supplied", {
  coords <- make_points(8)
  gd     <- make_gd(8)
  out    <- suppressWarnings(
    gdist.prep(n.Pops = 8, response = gd, samples = coords)
  )
  expect_true(all(c("gd", "pop") %in% names(out$df)))
  expect_equal(nrow(out$df), choose(8, 2))
})

test_that("gdist.prep attaches MLPE dyad metadata that matches legacy ZZ fitting", {
  coords <- make_points(6)
  gd <- make_gd(6)
  out <- suppressWarnings(
    gdist.prep(n.Pops = 6, response = gd, samples = coords, method = "costDistance")
  )

  expect_named(attr(out$df, "mlpe_pairs"), "pop")

  dat_attr <- out$df
  dat_attr$cd <- scale(seq_len(nrow(dat_attr)))

  fit_attr <- ResistanceGA2::mlpe_rga(
    formula = out$formula,
    data = dat_attr,
    ZZ = out$ZZ,
    REML = FALSE
  )

  dat_legacy <- dat_attr
  attr(dat_legacy, "mlpe_pairs") <- NULL

  fit_legacy <- ResistanceGA2::mlpe_rga(
    formula = out$formula,
    data = dat_legacy,
    ZZ = out$ZZ,
    REML = FALSE
  )

  expect_equal(as.numeric(logLik(fit_attr)),
               as.numeric(logLik(fit_legacy)),
               tolerance = 1e-8)
  expect_equal(unname(lme4::fixef(fit_attr)),
               unname(lme4::fixef(fit_legacy)),
               tolerance = 1e-8)
})

test_that("gdist.prep df is NULL when no response supplied", {
  coords <- make_points(8)
  out    <- gdist.prep(n.Pops = 8, samples = coords)
  expect_null(out$df)
})

# ----- Input validation errors ------------------------------------------------

test_that("gdist.prep errors when n.Pops != nrow(samples)", {
  coords <- make_points(10)
  expect_error(gdist.prep(n.Pops = 5, samples = coords), "n.Pops")
})

test_that("gdist.prep errors when response is a matrix", {
  coords <- make_points(10)
  expect_error(
    gdist.prep(n.Pops = 10, response = matrix(1:45, 9), samples = coords),
    "vector"
  )
})

test_that("gdist.prep errors when covariates is not a data frame", {
  coords <- make_points(10)
  gd     <- make_gd(10)
  expect_error(
    gdist.prep(n.Pops = 10, response = gd, samples = coords,
               covariates = 1:45),
    "data frame"
  )
})

test_that("gdist.prep errors on mismatched response and covariates lengths", {
  coords <- make_points(10)
  gd     <- make_gd(10)
  cov    <- data.frame(x = runif(10))  # wrong length: need choose(10,2) = 45
  expect_error(
    gdist.prep(n.Pops = 10, response = gd, samples = coords, covariates = cov),
    "same number"
  )
})

test_that("gdist.prep errors on unsupported min.max_dist", {
  coords <- make_points(10)
  gd     <- make_gd(10)
  expect_error(
    suppressWarnings(
      gdist.prep(n.Pops = 10, response = gd, samples = coords,
                 min.max_dist = c(0, 50))
    ),
    "not yet supported"
  )
})

# ----- Input format flexibility -----------------------------------------------

test_that("gdist.prep accepts SpatVector samples", {
  coords <- make_coords(10)
  sv     <- terra::vect(coords, type = "points")
  gd     <- make_gd(10)
  out    <- suppressWarnings(
    gdist.prep(n.Pops = 10, response = gd, samples = sv)
  )
  expect_equal(as.numeric(out$samples), as.numeric(coords))
})

test_that("gdist.prep errors on matrix samples", {
  coords <- make_coords(10)
  gd     <- make_gd(10)
  expect_error(
    gdist.prep(n.Pops = 10, response = gd, samples = coords),
    "SpatVector"
  )
})

test_that("gdist.prep errors on data.frame samples", {
  coords <- as.data.frame(make_coords(10))
  gd <- make_gd(5)
  expect_error(
    gdist.prep(n.Pops = 10, response = make_gd(10), samples = coords),
    "SpatVector"
  )
})

test_that("gdist.prep errors on invalid samples type", {
  gd <- make_gd(5)
  expect_error(
    gdist.prep(n.Pops = 5, response = gd, samples = list(x = 1:5, y = 1:5)),
    "SpatVector"
  )
})

# ----- Covariates and formula -------------------------------------------------

test_that("gdist.prep incorporates covariates into df", {
  coords <- make_points(8)
  gd     <- make_gd(8)
  n_pairs <- choose(8, 2)
  cov    <- data.frame(elevation = runif(n_pairs))
  out    <- suppressWarnings(
    gdist.prep(n.Pops = 8, response = gd, samples = coords,
               covariates = cov)
  )
  expect_true("elevation" %in% names(out$df))
})

test_that("gdist.prep updates formula when covariates and formula supplied", {
  coords <- make_points(8)
  gd     <- make_gd(8)
  n_pairs <- choose(8, 2)
  cov    <- data.frame(elev = runif(n_pairs))
  fml    <- as.formula("response ~ elev")
  out    <- suppressWarnings(
    gdist.prep(n.Pops = 8, response = gd, samples = coords,
               covariates = cov, formula = fml)
  )
  # Updated formula should contain cd
  expect_true(grepl("cd", deparse(out$formula)))
})

# ----- Warning on inverted genetic distance -----------------------------------

test_that("gdist.prep warns when gd decreases with Euclidean distance", {
  n      <- 8
  coords <- matrix(c(seq_len(n), rep(0, n)), ncol = 2)
  pts    <- terra::vect(coords, type = "points")
  ed     <- as.vector(dist(coords))
  # Invert: large Euclidean → small genetic distance
  gd     <- max(ed) - ed + rnorm(length(ed), 0, 0.01)
  expect_warning(
    gdist.prep(n.Pops = n, response = gd, samples = pts),
    "decreases"
  )
})
