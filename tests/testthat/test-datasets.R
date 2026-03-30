test_that("precomputed vignette outputs load with expected structure", {
  data_env <- new.env(parent = emptyenv())
  utils::data(
    "vignette_pairwise_outputs",
    package = "ResistanceGA2",
    envir = data_env
  )
  utils::data(
    "vignette_current_map",
    package = "ResistanceGA2",
    envir = data_env
  )

  vignette_pairwise_outputs <- get("vignette_pairwise_outputs", envir = data_env)
  current_map <- get("vignette_current_map", envir = data_env)

  expect_named(
    vignette_pairwise_outputs,
    c("gdistance_commute", "circuitscape_matrix")
  )

  expect_true(is.matrix(vignette_pairwise_outputs$gdistance_commute))
  expect_equal(dim(vignette_pairwise_outputs$gdistance_commute), c(25L, 25L))
  expect_true(is.numeric(vignette_pairwise_outputs$circuitscape_matrix))
  expect_equal(length(vignette_pairwise_outputs$circuitscape_matrix), choose(25L, 2L))

  if (inherits(current_map, "PackedSpatRaster")) {
    current_map <- terra::unwrap(current_map)
  }
  expect_true(inherits(current_map, "SpatRaster"))
  expect_equal(terra::nlyr(current_map), 1L)
  expect_identical(names(current_map), "continuous")
})
