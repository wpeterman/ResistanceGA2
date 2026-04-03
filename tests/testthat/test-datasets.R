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

test_that("rga2 demo data load with expected structure", {
  data_env <- new.env(parent = emptyenv())
  utils::data("rga2_demo", package = "ResistanceGA2", envir = data_env)
  utils::data("rga2_demo_covariates", package = "ResistanceGA2", envir = data_env)
  utils::data("rga2_demo_genetic", package = "ResistanceGA2", envir = data_env)

  rga2_demo <- get("rga2_demo", envir = data_env)
  rga2_demo_covariates <- get("rga2_demo_covariates", envir = data_env)
  rga2_demo_genetic <- get("rga2_demo_genetic", envir = data_env)

  expect_named(
    rga2_demo,
    c("sites", "sample_coords", "site_environment", "genetic_matrix", "genetic_vector")
  )

  expect_s3_class(rga2_demo$sites, "data.frame")
  expect_equal(nrow(rga2_demo$sites), 12L)
  expect_true(all(c("demo_id", "SiteName", "x", "y") %in% names(rga2_demo$sites)))

  expect_true(is.matrix(rga2_demo$sample_coords))
  expect_equal(dim(rga2_demo$sample_coords), c(12L, 2L))
  expect_identical(colnames(rga2_demo$sample_coords), c("x", "y"))

  expect_s3_class(rga2_demo$site_environment, "data.frame")
  expect_equal(nrow(rga2_demo$site_environment), 12L)
  expect_true(all(c("demo_id", "cti", "err27", "ffp", "gsp", "hli", "nlcd") %in%
    names(rga2_demo$site_environment)))

  expect_true(is.matrix(rga2_demo_genetic))
  expect_equal(dim(rga2_demo_genetic), c(12L, 12L))
  expect_identical(rownames(rga2_demo_genetic), colnames(rga2_demo_genetic))
  expect_true(isSymmetric(rga2_demo_genetic))

  expect_identical(rga2_demo$genetic_matrix, rga2_demo_genetic)
  expect_identical(rownames(rga2_demo$sample_coords), colnames(rga2_demo_genetic))
  expect_identical(rga2_demo$site_environment$demo_id, colnames(rga2_demo_genetic))
  expect_equal(length(rga2_demo$genetic_vector), choose(12L, 2L))
  expect_equal(rga2_demo$genetic_vector, lower(rga2_demo_genetic))

  if (inherits(rga2_demo_covariates, "PackedSpatRaster")) {
    rga2_demo_covariates <- terra::unwrap(rga2_demo_covariates)
  }
  expect_true(inherits(rga2_demo_covariates, "SpatRaster"))
  expect_equal(terra::nlyr(rga2_demo_covariates), 6L)
  expect_identical(
    names(rga2_demo_covariates),
    c("cti", "err27", "ffp", "gsp", "hli", "nlcd")
  )
})
