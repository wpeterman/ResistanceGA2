make_optimizer_fixture <- function() {
  cont <- terra::rast(
    nrows = 20,
    ncols = 20,
    xmin = 0,
    xmax = 10,
    ymin = 0,
    ymax = 10
  )
  set.seed(42)
  terra::values(cont) <- runif(terra::ncell(cont), 1, 10)
  names(cont) <- "continuous"

  breaks <- stats::quantile(terra::values(cont), probs = c(1 / 3, 2 / 3), na.rm = TRUE)
  cat_r <- terra::classify(
    cont,
    rbind(
      c(-Inf, breaks[1], 1),
      c(breaks[1], breaks[2], 2),
      c(breaks[2], Inf, 3)
    )
  )
  names(cat_r) <- "categorical"

  coords <- matrix(
    c(1, 1,
      2, 8,
      4, 3,
      6, 9,
      8, 2,
      9, 8),
    ncol = 2,
    byrow = TRUE
  )
  pts <- terra::vect(coords, type = "points")
  keep <- rep(1L, choose(nrow(coords), 2))

  base_gdist <- ResistanceGA2::gdist.prep(
    n.Pops = nrow(coords),
    samples = pts,
    keep = keep,
    method = "costDistance"
  )

  list(
    continuous = cont,
    categorical = cat_r,
    multi = c(cat_r, cont),
    points = pts,
    keep = keep,
    base_gdist = base_gdist,
    n_pops = nrow(coords)
  )
}

make_optimizer_inputs <- function(response,
                                  fixture,
                                  covariates = NULL,
                                  formula = NULL) {
  ResistanceGA2::gdist.prep(
    n.Pops = fixture$n_pops,
    samples = fixture$points,
    response = response,
    keep = fixture$keep,
    method = "costDistance",
    covariates = covariates,
    formula = formula
  )
}

test_that("Resistance.Opt_single returns finite objective for continuous surfaces", {
  fixture <- make_optimizer_fixture()
  ga_inputs <- suppressMessages(
    ResistanceGA2::GA.prep(
      raster = fixture$continuous,
      Results.dir = tempfile(pattern = "opt-single-cont-"),
      monitor = FALSE,
      quiet = TRUE
    )
  )

  parm <- c(3, 2.5, 50)
  target_surface <- ResistanceGA2::Resistance.tran(
    transformation = "Monomolecular",
    shape = parm[2],
    max = parm[3],
    r = fixture$continuous
  )
  response <- ResistanceGA2::Run_gdistance(fixture$base_gdist, target_surface)
  gdist_inputs <- make_optimizer_inputs(response, fixture)

  obj <- ResistanceGA2:::Resistance.Opt_single(
    PARM = parm,
    Resistance = fixture$continuous,
    gdist.inputs = gdist_inputs,
    GA.inputs = ga_inputs,
    iter = 1,
    quiet = TRUE
  )

  expect_true(is.finite(obj))
  expect_gt(obj, -99999)
})

test_that("Resistance.Opt_single rejects negative fitted relationships", {
  fixture <- make_optimizer_fixture()
  ga_inputs <- suppressMessages(
    ResistanceGA2::GA.prep(
      raster = fixture$continuous,
      Results.dir = tempfile(pattern = "opt-single-neg-"),
      monitor = FALSE,
      quiet = TRUE,
      select.trans = list("A")
    )
  )

  parm <- c(3, 2.5, 50)
  target_surface <- ResistanceGA2::Resistance.tran(
    transformation = "Monomolecular",
    shape = parm[2],
    max = parm[3],
    r = fixture$continuous
  )
  response <- ResistanceGA2::Run_gdistance(fixture$base_gdist, target_surface)

  expect_warning(
    gdist_inputs <- make_optimizer_inputs(-response, fixture),
    "decreases"
  )

  obj <- ResistanceGA2:::Resistance.Opt_single(
    PARM = parm,
    Resistance = fixture$continuous,
    gdist.inputs = gdist_inputs,
    GA.inputs = ga_inputs,
    iter = 1,
    quiet = TRUE
  )

  expect_identical(obj, -99999)
})

test_that("Resistance.Opt_single rejects oversized Ricker shapes", {
  fixture <- make_optimizer_fixture()
  ga_inputs <- suppressMessages(
    ResistanceGA2::GA.prep(
      raster = fixture$continuous,
      Results.dir = tempfile(pattern = "opt-single-ricker-"),
      monitor = FALSE,
      quiet = TRUE,
      select.trans = list("A")
    )
  )

  parm_ok <- c(3, 2.5, 50)
  target_surface <- ResistanceGA2::Resistance.tran(
    transformation = "Monomolecular",
    shape = parm_ok[2],
    max = parm_ok[3],
    r = fixture$continuous
  )
  response <- ResistanceGA2::Run_gdistance(fixture$base_gdist, target_surface)
  gdist_inputs <- make_optimizer_inputs(response, fixture)

  obj <- ResistanceGA2:::Resistance.Opt_single(
    PARM = c(4, 6, 50),
    Resistance = fixture$continuous,
    gdist.inputs = gdist_inputs,
    GA.inputs = ga_inputs,
    iter = 1,
    quiet = TRUE
  )

  expect_identical(obj, -99999)
})

test_that("Resistance.Opt_single returns finite objective for categorical surfaces", {
  fixture <- make_optimizer_fixture()
  ga_inputs <- suppressMessages(
    ResistanceGA2::GA.prep(
      raster = fixture$categorical,
      Results.dir = tempfile(pattern = "opt-single-cat-"),
      monitor = FALSE,
      quiet = TRUE
    )
  )

  parm <- c(1, 5, 10)
  df <- data.frame(id = terra::unique(fixture$categorical)[[1]], parm = parm)
  target_surface <- terra::subst(fixture$categorical, from = df$id, to = df$parm)
  response <- ResistanceGA2::Run_gdistance(fixture$base_gdist, target_surface)
  gdist_inputs <- make_optimizer_inputs(response, fixture)

  obj <- ResistanceGA2:::Resistance.Opt_single(
    PARM = parm,
    Resistance = fixture$categorical,
    gdist.inputs = gdist_inputs,
    GA.inputs = ga_inputs,
    iter = 1,
    quiet = TRUE
  )

  expect_true(is.finite(obj))
  expect_gt(obj, -99999)
})

test_that("Resistance.Opt_single.cov returns finite objective with covariates", {
  fixture <- make_optimizer_fixture()
  ga_inputs <- suppressMessages(
    ResistanceGA2::GA.prep(
      raster = fixture$continuous,
      Results.dir = tempfile(pattern = "opt-single-cov-"),
      monitor = FALSE,
      quiet = TRUE
    )
  )

  parm <- c(3, 2.5, 50)
  target_surface <- ResistanceGA2::Resistance.tran(
    transformation = "Monomolecular",
    shape = parm[2],
    max = parm[3],
    r = fixture$continuous
  )
  response_base <- ResistanceGA2::Run_gdistance(fixture$base_gdist, target_surface)
  covariates <- data.frame(elev = seq_along(response_base) / 10)
  response <- response_base + covariates$elev

  gdist_inputs <- make_optimizer_inputs(
    response = response,
    fixture = fixture,
    covariates = covariates,
    formula = response ~ elev
  )

  obj <- ResistanceGA2:::Resistance.Opt_single.cov(
    PARM = parm,
    Resistance = fixture$continuous,
    gdist.inputs = gdist_inputs,
    GA.inputs = ga_inputs,
    iter = 1,
    quiet = TRUE
  )

  expect_true(is.finite(obj))
  expect_gt(obj, -99999)
})

test_that("Resistance.Opt_multi returns finite objective for multisurface models", {
  fixture <- make_optimizer_fixture()
  ga_inputs <- suppressMessages(
    ResistanceGA2::GA.prep(
      raster = fixture$multi,
      Results.dir = tempfile(pattern = "opt-multi-"),
      monitor = FALSE,
      quiet = TRUE
    )
  )

  parm <- c(1, 5, 10, 3, 2.5, 50)
  target_surface <- ResistanceGA2::Combine_Surfaces(
    PARM = parm,
    GA.inputs = ga_inputs,
    rescale = FALSE
  )
  response <- ResistanceGA2::Run_gdistance(fixture$base_gdist, target_surface)
  gdist_inputs <- make_optimizer_inputs(response, fixture)

  obj <- ResistanceGA2:::Resistance.Opt_multi(
    PARM = parm,
    gdist.inputs = gdist_inputs,
    GA.inputs = ga_inputs,
    quiet = TRUE
  )

  expect_true(is.finite(obj))
  expect_gt(obj, -99999)
})

test_that("Resistance.Opt_multi.cov returns finite objective for multisurface covariate models", {
  fixture <- make_optimizer_fixture()
  ga_inputs <- suppressMessages(
    ResistanceGA2::GA.prep(
      raster = fixture$multi,
      Results.dir = tempfile(pattern = "opt-multi-cov-"),
      monitor = FALSE,
      quiet = TRUE
    )
  )

  parm <- c(1, 5, 10, 3, 2.5, 50)
  target_surface <- ResistanceGA2::Combine_Surfaces(
    PARM = parm,
    GA.inputs = ga_inputs,
    rescale = FALSE
  )
  response_base <- ResistanceGA2::Run_gdistance(fixture$base_gdist, target_surface)
  covariates <- data.frame(elev = seq_along(response_base) / 10)
  response <- response_base + covariates$elev

  gdist_inputs <- make_optimizer_inputs(
    response = response,
    fixture = fixture,
    covariates = covariates,
    formula = response ~ elev
  )

  obj <- ResistanceGA2:::Resistance.Opt_multi.cov(
    PARM = parm,
    gdist.inputs = gdist_inputs,
    GA.inputs = ga_inputs,
    quiet = TRUE
  )

  expect_true(is.finite(obj))
  expect_gt(obj, -99999)
})
