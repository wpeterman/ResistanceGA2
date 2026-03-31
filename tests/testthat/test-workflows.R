get_pkg_data <- function(name) {
  data_env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "ResistanceGA2", envir = data_env)
  get(name, envir = data_env)
}

unwrap_raster <- function(x) {
  if (inherits(x, "PackedSpatRaster")) {
    terra::unwrap(x)
  } else {
    x
  }
}

make_temp_dir <- function(prefix = "ResistanceGA2-test-") {
  path <- tempfile(pattern = prefix)
  dir.create(path, recursive = TRUE)
  paste0(normalizePath(path, winslash = "/"), "/")
}

make_pairwise_example <- function(n = 5, seed = 1) {
  set.seed(seed)
  coords <- matrix(runif(n * 2), ncol = 2)
  genetic_mat <- as.matrix(dist(coords))
  resistance_mat <- as.matrix(dist((coords * 1.25) + 0.1))

  list(
    coords = coords,
    genetic_mat = genetic_mat,
    genetic_vec = ResistanceGA2::lower(genetic_mat),
    resistance_mat = resistance_mat,
    resistance_vec = ResistanceGA2::lower(resistance_mat)
  )
}

make_smoke_inputs <- function(raster, prefix) {
  samples_df <- get_pkg_data("samples")
  pts <- terra::vect(samples_df[, 2:3], type = "points")
  response <- ResistanceGA2::lower(as.matrix(dist(samples_df[, 2:3])))

  list(
    gdist = ResistanceGA2::gdist.prep(
      n.Pops = nrow(samples_df),
      response = response,
      samples = pts,
      method = "costDistance"
    ),
    ga = ResistanceGA2::GA.prep(
      raster = raster,
      Results.dir = tempfile(pattern = prefix),
      pop.size = 10,
      maxiter = 1,
      run = 1,
      seed = 1,
      monitor = FALSE,
      quiet = TRUE
    )
  )
}

make_covariate_smoke_inputs <- function(raster, prefix, n = 8) {
  samples_df <- get_pkg_data("samples")[seq_len(n), , drop = FALSE]
  pts <- terra::vect(samples_df[, 2:3], type = "points")
  response_base <- ResistanceGA2::lower(as.matrix(dist(samples_df[, 2:3])))
  covariates <- data.frame(elev = seq_along(response_base) / 100)
  response <- response_base + covariates$elev

  list(
    gdist = ResistanceGA2::gdist.prep(
      n.Pops = nrow(samples_df),
      response = response,
      samples = pts,
      method = "costDistance",
      covariates = covariates,
      formula = response ~ elev
    ),
    ga = ResistanceGA2::GA.prep(
      raster = raster,
      Results.dir = tempfile(pattern = prefix),
      pop.size = 10,
      maxiter = 1,
      run = 1,
      seed = 1,
      monitor = FALSE,
      quiet = TRUE
    )
  )
}

test_that("MLPE.lmm fits vector, matrix, dist, and file inputs", {
  example <- make_pairwise_example()

  fit_vec <- ResistanceGA2::MLPE.lmm(
    resistance = example$resistance_vec,
    pairwise.genetic = example$genetic_vec,
    REML = FALSE,
    scale = FALSE
  )

  fit_mat <- ResistanceGA2::MLPE.lmm(
    resistance = example$resistance_mat,
    pairwise.genetic = example$genetic_vec,
    REML = FALSE
  )

  fit_dist <- ResistanceGA2::MLPE.lmm(
    resistance = stats::as.dist(example$resistance_mat),
    pairwise.genetic = example$genetic_vec,
    REML = FALSE
  )

  resistance_file <- file.path(make_temp_dir("mlpe-file-"), "resistance.txt")
  write.table(
    example$resistance_mat,
    file = resistance_file,
    row.names = TRUE,
    col.names = NA
  )

  fit_file <- ResistanceGA2::MLPE.lmm(
    resistance = resistance_file,
    pairwise.genetic = example$genetic_vec,
    REML = FALSE
  )

  expect_s4_class(fit_vec, "lmerMod")
  expect_s4_class(fit_mat, "lmerMod")
  expect_s4_class(fit_dist, "lmerMod")
  expect_s4_class(fit_file, "lmerMod")
  expect_true(all(is.finite(c(
    AIC(fit_vec),
    AIC(fit_mat),
    AIC(fit_dist),
    AIC(fit_file)
  ))))
})

test_that("Plot.trans returns plot objects and can write TIFF output", {
  skip_if_not(capabilities("tiff"))

  resistance_surfaces <- unwrap_raster(get_pkg_data("resistance_surfaces"))

  out_dir <- make_temp_dir("plot-trans-")

  p_numeric <- ResistanceGA2::Plot.trans(
    PARM = c(2.5, 100),
    Resistance = c(0, 1),
    transformation = "Monomolecular",
    marginal.plot = FALSE
  )

  p_raster <- ResistanceGA2::Plot.trans(
    PARM = c(2.5, 100),
    Resistance = terra::subset(resistance_surfaces, "continuous"),
    transformation = "Ricker",
    marginal.plot = TRUE,
    print.dir = out_dir,
    Name = "plot-test"
  )

  expect_s3_class(p_numeric, "ggplot")
  expect_true(any(c("ggExtraPlot", "ggplot") %in% class(p_raster)))
  expect_true(file.exists(file.path(out_dir, "Ricker_Transformation_plot-test.tif")))
})

test_that("Plot.trans validates Resistance input", {
  expect_error(
    ResistanceGA2::Plot.trans(
      PARM = c(2.5, 100),
      Resistance = 1:3,
      transformation = "Monomolecular",
      marginal.plot = FALSE
    ),
    "SpatRaster or a numeric vector of length 2"
  )
})

test_that("Diagnostic.Plots writes output for vector inputs and accepts file paths", {
  skip_if_not(capabilities("tiff"))

  example <- make_pairwise_example()
  out_dir <- make_temp_dir("diag-plots-")

  diag_out <- ResistanceGA2::Diagnostic.Plots(
    resistance.mat = example$resistance_vec,
    genetic.dist = example$genetic_vec,
    plot.dir = out_dir,
    type = "continuous",
    name = "diag-test"
  )

  expect_type(diag_out, "list")
  expect_true(file.exists(file.path(out_dir, "diag-test_DiagnosticPlots.tif")))

  resistance_file <- file.path(out_dir, "diag-resistance.txt")
  write.table(
    example$resistance_mat,
    file = resistance_file,
    row.names = TRUE,
    col.names = NA
  )

  expect_null(
    ResistanceGA2::Diagnostic.Plots(
      resistance.mat = resistance_file,
      genetic.dist = example$genetic_vec,
      plot.dir = out_dir,
      type = "continuous",
      name = "diag-file"
    )
  )
})

test_that("Combine_Surfaces returns a raster and percent contribution summary", {
  resistance_surfaces <- unwrap_raster(get_pkg_data("resistance_surfaces"))

  ga_inputs <- ResistanceGA2::GA.prep(
    raster = terra::subset(resistance_surfaces, c("categorical", "continuous")),
    Results.dir = tempfile(pattern = "combine-surfaces-"),
    monitor = FALSE,
    quiet = TRUE
  )

  parm <- c(1, 5, 10, 3, 2.5, 100)

  contribution <- ResistanceGA2::Combine_Surfaces(
    PARM = parm,
    GA.inputs = ga_inputs,
    p.contribution = TRUE
  )

  expect_named(contribution, c("percent.contribution", "combined.surface"))
  expect_true(inherits(contribution$combined.surface, "SpatRaster"))
  expect_equal(nrow(contribution$percent.contribution), 2L)

  out_dir <- make_temp_dir("combine-surfaces-out-")
  combined <- ResistanceGA2::Combine_Surfaces(
    PARM = parm,
    GA.inputs = ga_inputs,
    out = out_dir,
    File.name = "combined",
    rescale = FALSE
  )

  expect_true(inherits(combined, "SpatRaster"))
  expect_true(file.exists(file.path(out_dir, "combined.asc")))
})

test_that("Resist.boot ranks simple candidate distance matrices", {
  example <- make_pairwise_example()
  alt_mat <- as.matrix(dist((example$coords * 0.8) + 0.25))

  boot_out <- ResistanceGA2::Resist.boot(
    mod.names = c("model_1", "model_2"),
    dist.mat = list(example$resistance_mat, alt_mat),
    n.parameters = c(2, 2),
    sample.prop = 0.8,
    iters = 2,
    obs = nrow(example$genetic_mat),
    rank.method = "LL",
    genetic.mat = example$genetic_mat
  )

  expect_s3_class(boot_out, "data.frame")
  expect_equal(nrow(boot_out), 2L)
  expect_true(all(c("surface", "avg.rank", "Percent.top", "k") %in% names(boot_out)))
})

test_that("SS_optim completes a minimal gdistance workflow", {
  resistance_surfaces <- unwrap_raster(get_pkg_data("resistance_surfaces"))
  inputs <- make_smoke_inputs(
    raster = terra::subset(resistance_surfaces, "continuous"),
    prefix = "ss-optim-smoke-"
  )

  ss_out <- ResistanceGA2::SS_optim(
    gdist.inputs = inputs$gdist,
    GA.inputs = inputs$ga,
    dist_mod = FALSE,
    null_mod = FALSE,
    diagnostic_plots = FALSE
  )

  expect_named(
    ss_out,
    c("ContinuousResults", "CategoricalResults", "AICc", "MLPE",
      "Run.Time", "MLPE.list", "cd", "k", "ga")
  )
  expect_s3_class(ss_out$AICc, "data.frame")
  expect_length(ss_out$ga, 1L)
  expect_true("continuous" %in% names(ss_out$cd))
})

test_that("MS_optim completes a minimal multisurface workflow", {
  resistance_surfaces <- unwrap_raster(get_pkg_data("resistance_surfaces"))
  inputs <- make_smoke_inputs(
    raster = terra::subset(resistance_surfaces, c("categorical", "continuous")),
    prefix = "ms-optim-smoke-"
  )

  ms_out <- ResistanceGA2::MS_optim(
    gdist.inputs = inputs$gdist,
    GA.inputs = inputs$ga,
    diagnostic_plots = FALSE
  )

  expect_named(
    ms_out,
    c("GA.summary", "MLPE.model", "MLPE.model_REML", "AICc.tab",
      "cd", "percent.contribution", "k")
  )
  expect_s4_class(ms_out$GA.summary, "ga")
  expect_s3_class(ms_out$AICc.tab, "data.frame")
  expect_s3_class(ms_out$percent.contribution, "data.frame")
})

test_that("SS_optim covers covariate, distance, and null-model branches", {
  resistance_surfaces <- unwrap_raster(get_pkg_data("resistance_surfaces"))
  inputs <- make_covariate_smoke_inputs(
    raster = terra::subset(resistance_surfaces, c("categorical", "continuous")),
    prefix = "ss-optim-cov-"
  )

  ss_out <- ResistanceGA2::SS_optim(
    gdist.inputs = inputs$gdist,
    GA.inputs = inputs$ga,
    dist_mod = TRUE,
    null_mod = TRUE,
    diagnostic_plots = FALSE
  )

  expect_s3_class(ss_out$AICc, "data.frame")
  expect_true(all(c("Distance", "Null") %in% ss_out$AICc$Surface))
  expect_true(all(c("categorical", "continuous") %in% ss_out$AICc$Surface))
})

test_that("MS_optim covers the covariate branch", {
  resistance_surfaces <- unwrap_raster(get_pkg_data("resistance_surfaces"))
  inputs <- make_covariate_smoke_inputs(
    raster = terra::subset(resistance_surfaces, c("categorical", "continuous")),
    prefix = "ms-optim-cov-"
  )

  ms_out <- ResistanceGA2::MS_optim(
    gdist.inputs = inputs$gdist,
    GA.inputs = inputs$ga,
    diagnostic_plots = FALSE
  )

  expect_s3_class(ms_out$AICc.tab, "data.frame")
  expect_true(all(c("AICc", "LL", "R2m") %in% names(ms_out$AICc.tab)))
  expect_true(nrow(ms_out$percent.contribution) >= 2L)
})
