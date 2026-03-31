# Helper: generate a small MLPE-ready data frame
make_mlpe_data <- function(n = 5, seed = 1) {
  set.seed(seed)
  m_y <- matrix(rnorm(n^2), n, n)
  m_x <- matrix(rnorm(n^2), n, n)
  id  <- ResistanceGA2::To.From.ID(n)
  data.frame(
    y   = ResistanceGA2::lower(m_y),
    x   = ResistanceGA2::lower(m_x),
    pop = id$pop1
  )
}

make_mlpe_fixture <- function() {
  coords <- matrix(
    c(0, 0,
      1, 0,
      0, 1,
      1, 1),
    ncol = 2,
    byrow = TRUE
  )
  resistance_mat <- as.matrix(dist(coords))
  resistance_vec <- ResistanceGA2::lower(resistance_mat)
  genetic <- resistance_vec + seq_along(resistance_vec) / 10

  file_path <- tempfile(fileext = ".txt")
  writeLines(
    c(
      "id 1 2 3 4",
      "1 0 1 1 1.41421356",
      "2 1 0 1.41421356 1",
      "3 1 1.41421356 0 1",
      "4 1.41421356 1 1 0"
    ),
    con = file_path
  )

  list(
    matrix = resistance_mat,
    vector = resistance_vec,
    dist = stats::as.dist(resistance_mat),
    genetic = genetic,
    file = file_path
  )
}

test_that("mlpe_rga returns a lmerMod object", {
  df  <- make_mlpe_data(5)
  fit <- ResistanceGA2::mlpe_rga(y ~ x + (1 | pop), data = df)
  expect_s4_class(fit, "lmerMod")
})

test_that("mlpe_rga AIC is finite when REML = FALSE", {
  df  <- make_mlpe_data(5)
  fit <- ResistanceGA2::mlpe_rga(y ~ x + (1 | pop), data = df, REML = FALSE)
  expect_true(is.finite(AIC(fit)))
})

test_that("mlpe_rga fixef has expected names", {
  df  <- make_mlpe_data(6)
  fit <- ResistanceGA2::mlpe_rga(y ~ x + (1 | pop), data = df)
  fe  <- lme4::fixef(fit)
  expect_named(fe, c("(Intercept)", "x"))
})

test_that("mlpe_rga REML and ML give different log-likelihoods", {
  df       <- make_mlpe_data(6)
  fit_ml   <- ResistanceGA2::mlpe_rga(y ~ x + (1 | pop), data = df, REML = FALSE)
  fit_reml <- ResistanceGA2::mlpe_rga(y ~ x + (1 | pop), data = df, REML = TRUE)
  expect_false(isTRUE(all.equal(
    as.numeric(logLik(fit_ml)),
    as.numeric(logLik(fit_reml))
  )))
})

test_that("mlpe_rga accepts formula as character string", {
  df <- make_mlpe_data(5)
  expect_no_error(ResistanceGA2::mlpe_rga("y ~ x + (1 | pop)", data = df))
})

test_that("mlpe_rga keep argument subsets observations", {
  df   <- make_mlpe_data(6)
  keep <- rep(c(1L, 0L), length.out = nrow(df))
  fit  <- ResistanceGA2::mlpe_rga(y ~ x + (1 | pop), data = df, keep = keep)
  expect_s4_class(fit, "lmerMod")
})

test_that("mlpe_rga with supplied ZZ matches auto-computed ZZ", {
  n   <- 5
  df  <- make_mlpe_data(n)
  id  <- ResistanceGA2::To.From.ID(n)
  ZZ  <- ResistanceGA2:::ZZ.mat(id)
  fit_auto <- ResistanceGA2::mlpe_rga(y ~ x + (1 | pop), data = df)
  fit_man  <- ResistanceGA2::mlpe_rga(y ~ x + (1 | pop), data = df, ZZ = ZZ)
  expect_equal(
    as.numeric(logLik(fit_auto)),
    as.numeric(logLik(fit_man))
  )
})

test_that("mlpe_rga Poisson glmer runs without error", {
  set.seed(42)
  n   <- 5
  id  <- ResistanceGA2::To.From.ID(n)
  df  <- data.frame(
    y   = rpois(choose(n, 2), lambda = 5),
    x   = rnorm(choose(n, 2)),
    pop = id$pop1
  )
  expect_no_error(
    ResistanceGA2::mlpe_rga(y ~ x + (1 | pop), data = df, family = poisson)
  )
})

test_that("MLPE.lmm accepts matrix, vector, dist, and file resistance inputs", {
  fixture <- make_mlpe_fixture()
  on.exit(unlink(fixture$file), add = TRUE)

  fit_mat <- ResistanceGA2::MLPE.lmm(fixture$matrix, fixture$genetic, REML = FALSE)
  fit_vec <- ResistanceGA2::MLPE.lmm(fixture$vector, fixture$genetic, REML = FALSE)
  fit_dist <- ResistanceGA2::MLPE.lmm(fixture$dist, fixture$genetic, REML = FALSE)
  fit_file <- ResistanceGA2::MLPE.lmm(fixture$file, fixture$genetic, REML = FALSE)

  expect_s4_class(fit_mat, "lmerMod")
  expect_s4_class(fit_vec, "lmerMod")
  expect_s4_class(fit_dist, "lmerMod")
  expect_s4_class(fit_file, "lmerMod")
})

test_that("MLPE.lmm respects scale = FALSE", {
  fixture <- make_mlpe_fixture()
  on.exit(unlink(fixture$file), add = TRUE)

  fit <- ResistanceGA2::MLPE.lmm(
    resistance = fixture$vector,
    pairwise.genetic = fixture$genetic,
    REML = FALSE,
    scale = FALSE
  )

  expect_equal(as.numeric(fit@frame$resistance), fixture$vector, tolerance = 1e-8)
})

test_that("MLPE.lmm2 accepts dist inputs with supplied ID and ZZ", {
  fixture <- make_mlpe_fixture()
  on.exit(unlink(fixture$file), add = TRUE)

  id <- ResistanceGA2::To.From.ID(4)
  zz <- ResistanceGA2:::ZZ.mat(id)
  fit <- ResistanceGA2:::MLPE.lmm2(
    resistance = fixture$dist,
    response = fixture$genetic,
    REML = FALSE,
    ID = id,
    ZZ = zz
  )

  expect_s4_class(fit, "lmerMod")
})

test_that("MLPE.lmm_coef reads distance-matrix outputs and returns coefficient tables", {
  fixture <- make_mlpe_fixture()
  coeff_dir <- tempfile("mlpe-coef-")
  dir.create(coeff_dir)
  file.create(file.path(coeff_dir, "surface_costDistance_distMat.csv"))
  on.exit(unlink(coeff_dir, recursive = TRUE, force = TRUE), add = TRUE)
  on.exit(unlink(fixture$file), add = TRUE)

  out <- suppressWarnings(
    testthat::with_mocked_bindings(
      ResistanceGA2:::MLPE.lmm_coef(
        resistance = coeff_dir,
        genetic.dist = fixture$genetic,
        method = "gd"
      ),
      read.csv = function(...) fixture$matrix,
      .package = "ResistanceGA2"
    )
  )

  expect_true(is.matrix(out) || is.data.frame(out))
  expect_true(any(rownames(out) %in% "surface"))
})
