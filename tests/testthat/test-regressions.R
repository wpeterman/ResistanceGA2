make_mock_ga_result <- function(solution, fitness_value = 1) {
  if (!methods::isClass("mock_ga_result")) {
    methods::setClass(
      "mock_ga_result",
      slots = c(solution = "matrix", fitnessValue = "numeric")
    )
  }

  methods::new(
    "mock_ga_result",
    solution = solution,
    fitnessValue = fitness_value
  )
}

make_surface_stack_regression <- function() {
  categorical <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 1,
                             ymin = 0, ymax = 1)
  terra::values(categorical) <- rep(c(1, 2, 3, 1), length.out = terra::ncell(categorical))

  continuous <- terra::rast(categorical)
  terra::values(continuous) <- seq_len(terra::ncell(continuous))

  out <- c(categorical, continuous)
  names(out) <- c("categorical", "continuous")
  out
}

test_that("To.From.ID supports spatial expanded IDs when pop_n is supplied", {
  sp_loc <- terra::vect(cbind(c(1, 5, 10), c(1, 5, 10)), type = "points")

  out <- suppressWarnings(
    To.From.ID(
      sampled_pops = 3,
      pop_n = c(1L, 2L, 1L),
      spLoc = sp_loc,
      nb = 6
    )
  )

  expect_named(
    out,
    c("pop1.ind", "pop2.ind", "pop1.pop", "pop2.pop", "cor.grp", "corr_", "pop")
  )
  expect_equal(nrow(out), 5L)
  expect_true(all(out$corr_ %in% c(0, 1)))
})

test_that("SS_optim handles categorical gdistance covariate island branch", {
  set.seed(42)

  categorical <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 1,
                             ymin = 0, ymax = 1)
  terra::values(categorical) <- rep(c(1, 2), length.out = terra::ncell(categorical))
  names(categorical) <- "categorical"

  coords <- matrix(
    c(0, 0,
      1, 0,
      0, 1,
      1, 1,
      0.5, 0.5),
    ncol = 2,
    byrow = TRUE
  )
  response <- seq_len(10) / 10 + rnorm(10, sd = 0.05)
  covariates <- data.frame(elev = seq_len(10))
  samples <- terra::vect(coords, type = "points")

  gdist_inputs <- suppressWarnings(
    gdist.prep(
      n.Pops = 5,
      response = response,
      samples = samples,
      covariates = covariates,
      formula = gd ~ elev
    )
  )

  results_dir <- tempfile("rga2-ss-optim-")
  dir.create(results_dir, recursive = TRUE)
  plots_dir <- file.path(results_dir, "Plots")
  dir.create(plots_dir, recursive = TRUE)
  on.exit(unlink(results_dir, recursive = TRUE, force = TRUE), add = TRUE)

  ga_inputs <- list(
    k.value = 2,
    n.layers = 1L,
    Resistance.stack = categorical,
    layer.names = "categorical",
    gaisl = TRUE,
    surface.type = "cat",
    population = NULL,
    selection = NULL,
    pcrossover = 0.85,
    pmutation = 0.1,
    crossover = NULL,
    Min.Max = "min",
    min.list = list(c(1, 1)),
    max.list = list(c(100, 100)),
    numIslands = 2L,
    migrationRate = 0.1,
    migrationInterval = 1L,
    optim = FALSE,
    optimArgs = NULL,
    parallel = FALSE,
    pop.size = 6L,
    maxiter = 1L,
    run = 1L,
    percent.elite = 1L,
    mutation = NULL,
    seed = 1L,
    monitor = FALSE,
    quiet = TRUE,
    Results.dir = paste0(
      normalizePath(results_dir, winslash = "/", mustWork = FALSE),
      "/"
    ),
    Plots.dir = paste0(
      normalizePath(plots_dir, winslash = "/", mustWork = FALSE),
      "/"
    ),
    max.cat = 100,
    method = "LL",
    parm.type = data.frame(type = "cat", n.parm = 2, name = "categorical")
  )

  out <- testthat::with_mocked_bindings(
    suppressMessages(
      suppressWarnings(
        SS_optim(
          gdist.inputs = gdist_inputs,
          GA.inputs = ga_inputs,
          dist_mod = FALSE,
          null_mod = FALSE,
          diagnostic_plots = FALSE
        )
      )
    ),
    gaisl = function(...) {
      make_mock_ga_result(matrix(c(1, 2), nrow = 1), fitness_value = 5)
    },
    Run_gdistance = function(...) seq_along(gdist_inputs$response),
    Diagnostic.Plots = function(...) NULL,
    writeRaster = function(...) NULL,
    write.table = function(...) NULL,
    saveRDS = function(...) NULL,
    MLPE.lmm_coef = function(...) data.frame(),
    .package = "ResistanceGA2"
  )

  expect_type(out, "list")
  expect_s3_class(out$CategoricalResults, "data.frame")
  expect_equal(out$CategoricalResults$Surface, "categorical")
  expect_equal(length(out$cd), 1L)
})

test_that("Run_gdistance converts SpatRaster internally for gdistance", {
  r <- terra::rast(nrows = 6, ncols = 6, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  terra::values(r) <- seq_len(terra::ncell(r))

  pts <- terra::vect(
    matrix(
      c(0.1, 0.1,
        0.9, 0.1,
        0.1, 0.9,
        0.9, 0.9),
      ncol = 2,
      byrow = TRUE
    ),
    type = "points"
  )

  inputs <- gdist.prep(n.Pops = 4, samples = pts, method = "costDistance")
  out <- Run_gdistance(inputs, r)

  expect_false(identical(out, -99999))
  expect_equal(dim(out), c(4L, 4L))
})

test_that("all_comb auto-creates and uses the requested results directory", {
  surfaces <- make_surface_stack_regression()

  ga_inputs <- list(
    Results.dir = NULL,
    Write.dir = NULL,
    Plots.dir = NULL,
    n.layers = 2L,
    layer.names = names(surfaces),
    pop.size = 4L,
    Resistance.stack = surfaces,
    gaisl = FALSE,
    inputs = list(
      min.cat = 1,
      max.cat = 100,
      max.cont = 1000,
      select.trans = list(NA, 1:8),
      method = "LL",
      k.value = 2,
      pop.mult = 15,
      percent.elite = 5,
      type = "real-valued",
      pcrossover = 0.85,
      pmutation = 0.125,
      maxiter = 2L,
      run = 1L,
      keepBest = TRUE,
      population = NULL,
      selection = NULL,
      crossover = NULL,
      mutation = NULL,
      pop.size = 4L,
      parallel = FALSE,
      optim = FALSE,
      optim.method = "L-BFGS-B",
      poptim = 0,
      pressel = 1,
      control = list(fnscale = -1, maxit = 10),
      hessian = FALSE,
      monitor = FALSE,
      seed = 1L,
      quiet = TRUE
    )
  )

  gdist_inputs <- list(
    n.Pops = 3L,
    response = c(1, 2, 3)
  )

  results_dir <- tempfile("rga2-all-comb-")
  on.exit(unlink(results_dir, recursive = TRUE, force = TRUE), add = TRUE)

  out <- testthat::with_mocked_bindings(
    suppressMessages(
      all_comb(
        gdist.inputs = gdist_inputs,
        GA.inputs = ga_inputs,
        results.dir = results_dir,
        max.combination = c(2, 2),
        iters = 2,
        replicate = 1
      )
    ),
    SS_optim = function(...) {
      stop("SS_optim should not be called when min.combination > 1.")
    },
    GA.prep = function(raster, ...) {
      if (!inherits(raster, "SpatRaster")) {
        stop("Expected terra::SpatRaster input.")
      }
      list()
    },
    MS_optim = function(...) {
      list(
        cd = list(c(1, 2, 3)),
        k = data.frame(surface = "categorical.continuous", k = 2),
        AICc.tab = data.frame(surface = "categorical.continuous", AICc = 1)
      )
    },
    Resist.boot = function(mod.names,
                           dist.mat,
                           n.parameters,
                           sample.prop,
                           iters,
                           obs,
                           genetic.mat,
                           keep = NULL) {
      data.frame(model = mod.names, support = 1)
    },
    .package = "ResistanceGA2"
  )

  expect_type(out, "list")
  expect_null(out$ss.results)
  expect_equal(out$all.k$surface, "categorical.continuous")
  expect_equal(names(out$all.cd), "categorical.continuous")
  expect_true(dir.exists(results_dir))
  expect_true(dir.exists(file.path(results_dir, "rep_1")))
})
