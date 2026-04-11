make_run_gdistance_fixture <- function(method = "costDistance", keep = NULL, ...) {
  r <- terra::rast(nrows = 6, ncols = 6, xmin = 0, xmax = 6, ymin = 0, ymax = 6)
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
      method = method,
      keep = keep,
      ...
    )
  )
}

make_reference_transition <- function(fixture, correction) {
  r_gd <- raster::raster(fixture$r)
  tr <- gdistance::transition(
    x = r_gd,
    transitionFunction = fixture$inputs$transitionFunction,
    directions = fixture$inputs$directions
  )
  gdistance::geoCorrection(tr, correction, scl = TRUE)
}

selected_pair_reference <- function(fixture, transition, distance_fun) {
  keep <- which(fixture$inputs$keep == 1)
  id <- as.matrix(fixture$inputs$ID)
  id <- apply(id, 2, function(x) as.integer(as.character(x)))

  vapply(
    keep,
    function(i) {
      pts <- as.integer(id[i, 1:2])
      as.numeric(distance_fun(transition, fixture$inputs$samples[pts, , drop = FALSE]))
    },
    numeric(1)
  )
}

test_that("Run_gdistance preserves gdistance costDistance output", {
  fixture <- make_run_gdistance_fixture(method = "costDistance")
  trC <- make_reference_transition(fixture, "c")

  expected <- gdistance::costDistance(trC, fixture$inputs$samples)
  observed <- ResistanceGA2::Run_gdistance(fixture$inputs, fixture$r)

  expect_s3_class(observed, "dist")
  expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 1e-8)
})

test_that("Run_gdistance preserves selected-pair costDistance output", {
  keep <- rep(0L, choose(5, 2))
  keep[c(1, 4, 8)] <- 1L
  fixture <- make_run_gdistance_fixture(method = "costDistance", keep = keep)
  trC <- make_reference_transition(fixture, "c")

  expected <- selected_pair_reference(fixture, trC, gdistance::costDistance)
  observed <- ResistanceGA2::Run_gdistance(fixture$inputs, fixture$r)

  expect_type(observed, "double")
  expect_equal(observed, expected, tolerance = 1e-8)
})

test_that("Run_gdistance preserves gdistance commuteDistance output", {
  fixture <- make_run_gdistance_fixture(method = "commuteDistance")
  trR <- make_reference_transition(fixture, "r")

  expected <- gdistance::commuteDistance(trR, fixture$inputs$samples) / 1000
  observed <- ResistanceGA2::Run_gdistance(fixture$inputs, fixture$r)

  expect_s3_class(observed, "dist")
  expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 1e-8)
})

test_that("Run_gdistance preserves selected-pair commuteDistance output", {
  keep <- rep(0L, choose(5, 2))
  keep[c(2, 5, 9)] <- 1L
  fixture <- make_run_gdistance_fixture(method = "commuteDistance", keep = keep)
  trR <- make_reference_transition(fixture, "r")

  expected <- selected_pair_reference(
    fixture,
    trR,
    function(x, coords) gdistance::commuteDistance(x, coords) / 1000
  )
  observed <- ResistanceGA2::Run_gdistance(fixture$inputs, fixture$r)

  expect_type(observed, "double")
  expect_equal(observed, expected, tolerance = 1e-8)
})

test_that("Run_gdistance preserves custom transitionFunction output", {
  fixture <- make_run_gdistance_fixture(method = "commuteDistance")
  fixture$inputs$transitionFunction <- function(x) 1 / max(x)
  trR <- make_reference_transition(fixture, "r")

  expected <- gdistance::commuteDistance(trR, fixture$inputs$samples) / 1000
  observed <- ResistanceGA2::Run_gdistance(fixture$inputs, fixture$r)

  expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 1e-8)
})

test_that("Run_gdistance preserves empty selected-pair output", {
  keep <- rep(0L, choose(5, 2))
  fixture <- make_run_gdistance_fixture(method = "costDistance", keep = keep)

  expect_identical(
    ResistanceGA2::Run_gdistance(fixture$inputs, fixture$r),
    numeric(0)
  )
})

test_that("Run_gdistance preserves commuteDistance output on a 250x250 raster", {
  r <- terra::rast(nrows = 250, ncols = 250, xmin = 0, xmax = 250,
                   ymin = 0, ymax = 250)
  set.seed(250)
  terra::values(r) <- runif(terra::ncell(r), 1, 100)

  coords <- cbind(
    runif(12, 1, 248),
    runif(12, 1, 248)
  )
  inputs <- ResistanceGA2::gdist.prep(
    n.Pops = nrow(coords),
    samples = terra::vect(coords, type = "points"),
    method = "commuteDistance"
  )

  tr <- gdistance::transition(
    raster::raster(r),
    inputs$transitionFunction,
    inputs$directions
  )
  trR <- gdistance::geoCorrection(tr, "r", scl = TRUE)

  expected <- gdistance::commuteDistance(trR, inputs$samples) / 1000
  observed <- ResistanceGA2::Run_gdistance(inputs, r)

  expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 1e-8)
})

test_that("Run_gdistance transition cache recomputes changed raster values", {
  fixture <- make_run_gdistance_fixture(method = "commuteDistance")
  observed_first <- ResistanceGA2::Run_gdistance(fixture$inputs, fixture$r)

  fixture$r2 <- fixture$r
  terra::values(fixture$r2) <- rev(terra::values(fixture$r2))
  trR <- make_reference_transition(
    list(r = fixture$r2, inputs = fixture$inputs),
    "r"
  )

  expected_second <- gdistance::commuteDistance(
    trR,
    fixture$inputs$samples
  ) / 1000
  observed_second <- ResistanceGA2::Run_gdistance(fixture$inputs, fixture$r2)

  expect_false(isTRUE(all.equal(as.numeric(observed_first), as.numeric(observed_second))))
  expect_equal(as.numeric(observed_second), as.numeric(expected_second), tolerance = 1e-8)
})

test_that("Run_gdistance aggregate approximation uses aggregated raster", {
  fixture <- make_run_gdistance_fixture(
    method = "commuteDistance",
    commute.approx = "aggregate",
    approx.factor = 2L,
    approx.scale = FALSE
  )
  exact_inputs <- fixture$inputs
  exact_inputs$commute.approx <- "none"
  r_agg <- terra::aggregate(fixture$r, fact = 2L, fun = mean, na.rm = TRUE)

  expected <- ResistanceGA2::Run_gdistance(exact_inputs, r_agg)
  observed <- ResistanceGA2::Run_gdistance(fixture$inputs, fixture$r)

  expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 1e-8)
})

test_that("Run_gdistance exact override ignores stored aggregate approximation", {
  fixture <- make_run_gdistance_fixture(
    method = "commuteDistance",
    commute.approx = "aggregate",
    approx.factor = 2L
  )
  exact_inputs <- fixture$inputs
  exact_inputs$commute.approx <- "none"

  expected <- ResistanceGA2::Run_gdistance(exact_inputs, fixture$r)
  observed <- ResistanceGA2::Run_gdistance(
    fixture$inputs,
    fixture$r,
    commute.approx = "none"
  )

  expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 1e-8)
})
