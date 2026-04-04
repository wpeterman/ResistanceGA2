test_that("site_pairwise_diff returns expected pairwise differences", {
  env <- data.frame(
    temp = c(10, 13, 15),
    precip = c(100, 80, 110)
  )

  out <- ResistanceGA2::site_pairwise_diff(env)

  expect_s3_class(out, "data.frame")
  expect_named(out, c("temp_diff", "precip_diff"))
  expect_equal(out$temp_diff, c(3, 5, 2))
  expect_equal(out$precip_diff, c(20, 10, 30))

  vec_out <- ResistanceGA2::site_pairwise_diff(c(2, 5, 9))
  expect_named(vec_out, "x_diff")
  expect_equal(vec_out$x_diff, c(3, 7, 4))
})

test_that("site_pairwise helpers support scaling and constant variables", {
  env <- data.frame(
    temp = c(1, 3, 5),
    constant = c(2, 2, 2)
  )

  diff_out <- ResistanceGA2::site_pairwise_diff(env, scale = TRUE)
  dist_out <- ResistanceGA2::site_pairwise_dist(env, scale = TRUE)

  expect_equal(diff_out$temp_diff, c(1, 2, 1), tolerance = 1e-8)
  expect_equal(diff_out$constant_diff, c(0, 0, 0), tolerance = 1e-8)
  expect_equal(dist_out$env_dist, c(1, 2, 1), tolerance = 1e-8)
})

test_that("site_pairwise helpers support keep and pop_n expansion", {
  env <- data.frame(
    temp = c(1, 5, 9),
    precip = c(4, 1, 7)
  )
  pop_n <- c(2L, 1L, 2L)

  expected_diff <- ResistanceGA2:::expand.mat(abs(outer(env$temp, env$temp, "-")), pop_n)
  expected_dist <- ResistanceGA2:::expand.mat(
    as.matrix(stats::dist(scale(env))),
    pop_n
  )

  diff_out <- ResistanceGA2::site_pairwise_diff(env["temp"], pop_n = pop_n)
  dist_out <- ResistanceGA2::site_pairwise_dist(env, pop_n = pop_n)

  expect_equal(diff_out$temp_diff, expected_diff, tolerance = 1e-8)
  expect_equal(dist_out$env_dist, expected_dist, tolerance = 1e-8)

  keep <- rep(c(1L, 0L), length.out = length(expected_diff))
  diff_keep <- ResistanceGA2::site_pairwise_diff(env["temp"], pop_n = pop_n, keep = keep)
  dist_keep <- ResistanceGA2::site_pairwise_dist(env, pop_n = pop_n, keep = keep)

  expect_equal(diff_keep$temp_diff, expected_diff[keep == 1], tolerance = 1e-8)
  expect_equal(dist_keep$env_dist, expected_dist[keep == 1], tolerance = 1e-8)
})

test_that("site_pairwise helpers validate site-level data", {
  expect_error(
    ResistanceGA2::site_pairwise_diff(data.frame(temp = c(1, NA, 3))),
    "Missing values"
  )

  expect_error(
    ResistanceGA2::site_pairwise_dist(data.frame(temp = c("a", "b", "c"))),
    "numeric"
  )

  expect_error(
    ResistanceGA2::site_pairwise_diff(data.frame(temp = 1:3), vars = "precip"),
    "vars"
  )

  expect_error(
    ResistanceGA2::site_pairwise_dist(data.frame(temp = 1:3), name = ""),
    "name"
  )
})

test_that("site_pairwise helpers integrate with gdist.prep covariates", {
  coords <- matrix(
    c(0, 0,
      1, 2,
      3, 1,
      4, 4,
      6, 2),
    ncol = 2,
    byrow = TRUE
  )
  pts <- terra::vect(coords, type = "points")
  response <- ResistanceGA2::lower(as.matrix(stats::dist(coords)))
  env <- data.frame(
    temp = c(10, 12, 15, 11, 13),
    precip = c(100, 120, 80, 110, 95)
  )

  covariates <- cbind(
    ResistanceGA2::site_pairwise_diff(env, vars = "temp"),
    ResistanceGA2::site_pairwise_dist(env, vars = c("temp", "precip"), name = "climate_dist")
  )

  out <- ResistanceGA2::gdist.prep(
    n.Pops = nrow(coords),
    response = response,
    samples = pts,
    covariates = covariates,
    formula = response ~ temp_diff + climate_dist,
    method = "costDistance"
  )

  expect_true(all(c("temp_diff", "climate_dist") %in% names(out$df)))
  expect_match(paste(deparse(out$formula), collapse = " "), "temp_diff")
  expect_match(paste(deparse(out$formula), collapse = " "), "climate_dist")
})
