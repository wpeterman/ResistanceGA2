# Helper: small 0-10 scaled SpatRaster
make_raster <- function(n = 25, seed = 1) {
  set.seed(seed)
  r <- terra::rast(nrows = n, ncols = n,
                   xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  terra::values(r) <- runif(n^2, 0, 10)
  r
}

# ----- Resistance.tran --------------------------------------------------------

test_that("Resistance.tran returns a SpatRaster", {
  r   <- make_raster()
  out <- Resistance.tran("Monomolecular", shape = 2, max = 100, r = r)
  expect_true(inherits(out, "SpatRaster"))
})

test_that("Resistance.tran by name == by number (Monomolecular = 3)", {
  r    <- make_raster()
  out3 <- Resistance.tran(3,                shape = 2, max = 100, r = r)
  outn <- Resistance.tran("Monomolecular",  shape = 2, max = 100, r = r)
  expect_equal(terra::values(out3), terra::values(outn))
})

test_that("Resistance.tran Distance (eq 9) sets all cells to 1", {
  r   <- make_raster()
  out <- Resistance.tran(9, shape = 1, max = 1, r = r)
  vals <- as.vector(terra::values(out))
  expect_true(all(vals == 1))
})

test_that("Resistance.tran all equations return finite surfaces", {
  r    <- make_raster()
  eqs  <- 1:8
  for (eq in eqs) {
    out  <- Resistance.tran(eq, shape = 2, max = 50, r = r)
    minv <- terra::global(out, "min", na.rm = TRUE)[[1]]
    maxv <- terra::global(out, "max", na.rm = TRUE)[[1]]
    expect_true(is.finite(minv), info = paste("eq =", eq, "min finite"))
    expect_true(is.finite(maxv), info = paste("eq =", eq, "max finite"))
    expect_true(maxv <= 1e6 + 1,  info = paste("eq =", eq, "max <= 1e6"))
  }
})

test_that("Resistance.tran Monomolecular values increase with input values", {
  r   <- make_raster()
  out <- Resistance.tran("Monomolecular", shape = 2, max = 100, r = r)
  # Monomolecular is monotonically increasing: higher input → higher output
  r_vals   <- as.vector(terra::values(r))
  out_vals <- as.vector(terra::values(out))
  expect_true(cor(r_vals, out_vals, use = "complete.obs") > 0)
})

test_that("Resistance.tran Reverse Monomolecular values decrease with input", {
  r   <- make_raster()
  out <- Resistance.tran("Reverse Monomolecular", shape = 2, max = 100, r = r)
  r_vals   <- as.vector(terra::values(r))
  out_vals <- as.vector(terra::values(out))
  expect_true(cor(r_vals, out_vals, use = "complete.obs") < 0)
})

# ----- SHNe -------------------------------------------------------------------

test_that("SHNe returns data frame with correct columns", {
  Ne  <- c(10, 50, 100, 25)
  out <- SHNe(n.samples = 4, pop.size = Ne)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("dhm", "di"))
  expect_equal(nrow(out), choose(4, 2))
})

test_that("SHNe output='matrix' returns list of two matrices", {
  Ne  <- c(10, 50, 100, 25)
  out <- SHNe(n.samples = 4, pop.size = Ne, output = "matrix")
  expect_type(out, "list")
  expect_named(out, c("dhm", "di"))
  expect_equal(dim(out$dhm), c(4L, 4L))
})

test_that("SHNe dhm and di are negative and positive respectively", {
  Ne  <- rep(c(10, 100), 3)
  out <- SHNe(n.samples = 6, pop.size = Ne)
  expect_true(all(out$dhm <= 0))
  expect_true(all(out$di  >= 0))
})
