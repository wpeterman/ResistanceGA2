make_gdistance_optimization_stack <- function() {
  categorical <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4,
                             ymin = 0, ymax = 4)
  terra::values(categorical) <- c(
    1, 1, 2, 2,
    1, 2, 2, 2,
    3, 3, 3, 4,
    3, 3, 4, 4
  )

  continuous <- terra::rast(categorical)
  terra::values(continuous) <- seq_len(terra::ncell(continuous))

  stack <- c(categorical, continuous)
  names(stack) <- c("categorical", "continuous")
  stack
}

make_gdistance_optimization_inputs <- function(stack) {
  list(
    Resistance.stack = stack,
    raster = stack,
    surface.type = c("cat", "cont")
  )
}

make_gdistance_optimization_run_fixture <- function(...) {
  r <- terra::rast(nrows = 6, ncols = 6, xmin = 0, xmax = 6,
                   ymin = 0, ymax = 6)
  terra::values(r) <- seq_len(terra::ncell(r)) + 5

  coords <- matrix(
    c(0.5, 5.5,
      1.5, 4.5,
      3.5, 3.5,
      4.5, 2.5,
      2.5, 1.5),
    ncol = 2,
    byrow = TRUE
  )

  list(
    r = r,
    inputs = ResistanceGA2::gdist.prep(
      n.Pops = nrow(coords),
      samples = terra::vect(coords, type = "points"),
      ...
    )
  )
}

test_that("gdistance optimization approximation aggregates layers by type", {
  stack <- make_gdistance_optimization_stack()
  ga_inputs <- make_gdistance_optimization_inputs(stack)
  gdist_inputs <- list(
    method = "commuteDistance",
    commute.approx = "aggregate",
    approx.factor = 2L,
    approx.scale = FALSE
  )

  out <- ResistanceGA2:::.rga_prepare_gdistance_optimization_inputs(
    ga_inputs,
    gdist_inputs
  )
  approx_stack <- out$gdistance.approx.Resistance.stack

  expect_s4_class(approx_stack, "SpatRaster")
  expect_equal(terra::nrow(approx_stack), 2L)
  expect_equal(terra::ncol(approx_stack), 2L)
  expect_identical(names(approx_stack), names(stack))
  expect_equal(terra::values(approx_stack[[1]])[, 1], c(1, 2, 3, 4))

  expected_continuous <- terra::aggregate(
    stack[[2]],
    fact = 2L,
    fun = mean,
    na.rm = TRUE
  )
  expect_equal(
    terra::values(approx_stack[[2]])[, 1],
    terra::values(expected_continuous)[, 1]
  )
  expect_identical(out$gdistance.approx.factor, 2L)
  expect_false(out$gdistance.approx.scale)
})

test_that("gdistance optimization approximation is opt-in for commuteDistance", {
  stack <- make_gdistance_optimization_stack()
  ga_inputs <- make_gdistance_optimization_inputs(stack)

  no_approx <- ResistanceGA2:::.rga_prepare_gdistance_optimization_inputs(
    ga_inputs,
    list(method = "commuteDistance", commute.approx = "none")
  )
  cost_approx <- ResistanceGA2:::.rga_prepare_gdistance_optimization_inputs(
    ga_inputs,
    list(
      method = "costDistance",
      commute.approx = "aggregate",
      approx.factor = 2L
    )
  )

  expect_null(no_approx$gdistance.approx.Resistance.stack)
  expect_null(cost_approx$gdistance.approx.Resistance.stack)
})

test_that("gdistance optimization runs on the prepared approximate raster", {
  fixture <- make_gdistance_optimization_run_fixture(
    method = "commuteDistance",
    commute.approx = "aggregate",
    approx.factor = 2L,
    approx.scale = TRUE
  )
  ga_inputs <- list(
    Resistance.stack = fixture$r,
    raster = fixture$r,
    surface.type = "cont"
  )
  ga_inputs <- ResistanceGA2:::.rga_prepare_gdistance_optimization_inputs(
    ga_inputs,
    fixture$inputs
  )
  r_approx <- ga_inputs$gdistance.approx.Resistance.stack

  exact_inputs <- fixture$inputs
  exact_inputs$commute.approx <- "none"
  expected <- ResistanceGA2::Run_gdistance(exact_inputs, r_approx) * 4
  observed <- ResistanceGA2:::.rga_run_gdistance_optimization(
    fixture$inputs,
    r_approx,
    ga_inputs
  )

  expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 1e-8)
})

test_that("gdistance exact helper ignores stored optimization approximation", {
  fixture <- make_gdistance_optimization_run_fixture(
    method = "commuteDistance",
    commute.approx = "aggregate",
    approx.factor = 2L
  )

  exact_inputs <- fixture$inputs
  exact_inputs$commute.approx <- "none"
  expected <- ResistanceGA2::Run_gdistance(exact_inputs, fixture$r)
  observed <- ResistanceGA2:::.rga_run_gdistance_exact(
    fixture$inputs,
    fixture$r
  )

  expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 1e-8)
})
