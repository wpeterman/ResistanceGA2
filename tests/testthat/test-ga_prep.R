make_surface_stack <- function() {
  categorical <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 1,
                             ymin = 0, ymax = 1)
  terra::values(categorical) <- rep(c(1, 2, 3, 1), length.out = terra::ncell(categorical))

  continuous <- terra::rast(categorical)
  terra::values(continuous) <- seq_len(terra::ncell(continuous))

  out <- c(categorical, continuous)
  names(out) <- c("categorical", "continuous")
  out
}

test_that("GA.prep uses raster inputs and normalizes output directories", {
  surfaces <- make_surface_stack()
  results_dir <- tempfile("rga2-results-")
  on.exit(unlink(results_dir, recursive = TRUE, force = TRUE), add = TRUE)

  out <- suppressMessages(
    GA.prep(
      raster = surfaces,
      Results.dir = results_dir,
      monitor = FALSE,
      quiet = TRUE
    )
  )

  expected_root <- paste0(
    sub("[/\\\\]+$", "", normalizePath(results_dir, winslash = "/", mustWork = FALSE)),
    "/"
  )

  expect_equal(out$n.layers, 2L)
  expect_identical(out$layer.names, c("categorical", "continuous"))
  expect_identical(out$surface.type, c("cat", "cont"))
  expect_s4_class(out$raster, "SpatRaster")
  expect_true(dir.exists(file.path(results_dir, "Results")))
  expect_true(dir.exists(file.path(results_dir, "Results", "Plots")))
  expect_identical(out$Results.dir, paste0(expected_root, "Results/"))
  expect_identical(out$Plots.dir, paste0(expected_root, "Results/Plots/"))
})

test_that("GA.prep requires select.trans to be provided as a list", {
  surfaces <- make_surface_stack()[[2]]
  results_dir <- tempfile("rga2-results-")
  on.exit(unlink(results_dir, recursive = TRUE, force = TRUE), add = TRUE)

  expect_error(
    GA.prep(
      raster = surfaces,
      Results.dir = results_dir,
      select.trans = "M",
      monitor = FALSE,
      quiet = TRUE
    ),
    "list"
  )
})

test_that("GA.prep rejects non-raster inputs", {
  results_dir <- tempfile("rga2-results-")
  on.exit(unlink(results_dir, recursive = TRUE, force = TRUE), add = TRUE)

  expect_error(
    GA.prep(
      raster = tempdir(),
      Results.dir = results_dir,
      monitor = FALSE,
      quiet = TRUE
    ),
    "SpatRaster"
  )
})

test_that("GA.prep configures all_comb mode, island defaults, and optimArgs", {
  surfaces <- make_surface_stack()

  out <- GA.prep(
    raster = surfaces,
    Results.dir = "all_comb",
    gaisl = TRUE,
    pop.size = NULL,
    run = NULL,
    numIslands = NULL,
    migrationRate = NULL,
    migrationInterval = NULL,
    optim = TRUE,
    monitor = FALSE,
    quiet = TRUE
  )

  expect_null(out$Results.dir)
  expect_null(out$Write.dir)
  expect_null(out$Plots.dir)
  expect_true(out$gaisl)
  expect_true(isTRUE(out$parallel))
  expect_equal(out$numIslands, 4)
  expect_equal(out$migrationRate, 0.10)
  expect_equal(out$migrationInterval, 10)
  expect_equal(out$run, 50)
  expect_equal(out$pop.size, 80)
  expect_type(out$optimArgs, "list")
  expect_equal(out$select.trans[[2]], c(1, 3, 5, 7, 9))
})

test_that("GA.prep warns and defaults Results.dir to the current working directory", {
  surfaces <- make_surface_stack()[[2]]
  original_wd <- getwd()
  test_wd <- tempfile("rga2-ga-prep-wd-")
  dir.create(test_wd, recursive = TRUE)
  on.exit({
    setwd(original_wd)
    unlink(test_wd, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  setwd(test_wd)

  out <- suppressMessages(
    expect_warning(
      GA.prep(
        raster = surfaces,
        Results.dir = NULL,
        monitor = FALSE,
        quiet = TRUE
      ),
      "not specified"
    )
  )

  root <- paste0(normalizePath(test_wd, winslash = "/", mustWork = TRUE), "/")
  expect_identical(out$Results.dir, paste0(root, "Results/"))
  expect_identical(out$Plots.dir, paste0(root, "Results/Plots/"))
})
