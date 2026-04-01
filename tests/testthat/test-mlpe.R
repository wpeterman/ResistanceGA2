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

legacy_mlpe_rga <- function(formula,
                            data,
                            REML = FALSE,
                            ZZ = NULL,
                            keep = NULL,
                            ...) {
  if (!inherits(formula, "formula")) {
    formula <- as.formula(formula)
  }

  if (is.null(ZZ)) {
    obs <- 0.5 * (sqrt((8 * nrow(data)) + 1) + 1)
    ID <- ResistanceGA2::To.From.ID(obs)
    ZZ <- ResistanceGA2:::ZZ.mat(ID)
  }

  if (!is.null(keep)) {
    data <- data[keep == 1, , drop = FALSE]
    ZZ <- ZZ[, keep == 1, drop = FALSE]
    pop <- ID$pop1[keep == 1]
    ZZ <- ZZ[rownames(ZZ) %in% unique(as.character(pop)), , drop = FALSE]
  }

  args <- list(...)

  if (any("family" %in% names(args))) {
    mod <- lme4::glFormula(formula, data = data, ...)
    mod$reTrms$Zt <- ZZ
    dfun <- do.call(lme4::mkGlmerDevfun, mod)
    opt <- lme4::optimizeGlmer(dfun)
  } else {
    mod <- lme4::lFormula(formula, data = data, REML = REML)
    mod$reTrms$Zt <- ZZ
    dfun <- do.call(lme4::mkLmerDevfun, mod)
    opt <- lme4::optimizeLmer(dfun)
  }

  lme4::mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr)
}

legacy_MLPE_lmm <- function(resistance,
                            pairwise.genetic,
                            REML = FALSE,
                            ID = NULL,
                            ZZ = NULL,
                            scale = TRUE) {
  response <- pairwise.genetic

  if (class(resistance)[[1]] == "dist") {
    mm <- as.vector(resistance)
    m <- attr(resistance, "Size")
    mm <- mm[which(mm != -1)]

    if (is.null(ID)) {
      ID <- ResistanceGA2::To.From.ID(m)
    }
    if (is.null(ZZ)) {
      ZZ <- ResistanceGA2:::ZZ.mat(ID = ID)
    }
  } else if (!is.character(resistance)) {
    if (is.vector(resistance)) {
      mm <- resistance
      m <- 0.5 * (sqrt((8 * length(mm)) + 1) + 1)
      mm <- mm[which(mm != -1)]
    } else {
      mm <- resistance
      m <- nrow(mm)
      mm <- ResistanceGA2::lower(mm)
      mm <- mm[which(mm != -1)]
    }

    if (is.null(ID)) {
      ID <- ResistanceGA2::To.From.ID(m)
    }
    if (is.null(ZZ)) {
      ZZ <- ResistanceGA2:::ZZ.mat(ID = ID)
    }
  } else {
    mm <- (read.table(resistance)[-1, -1])
    m <- nrow(mm)
    mm <- ResistanceGA2::lower(mm)
    mm <- mm[which(mm != -1)]

    if (is.null(ID)) {
      ID <- ResistanceGA2::To.From.ID(m)
    }
    if (is.null(ZZ)) {
      ZZ <- ResistanceGA2:::ZZ.mat(ID = ID)
    }
  }

  cs.matrix <- if (isTRUE(scale)) {
    base::scale(mm, center = TRUE, scale = TRUE)
  } else {
    mm
  }

  dat <- data.frame(ID, resistance = cs.matrix, response = response)
  colnames(dat) <- c("pop1", "pop2", "resistance", "response")

  mod <- lme4::lFormula(response ~ resistance + (1 | pop1),
                        data = dat,
                        REML = REML)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(lme4::mkLmerDevfun, mod)
  opt <- lme4::optimizeLmer(dfun)
  lme4::mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr)
}

summarize_mlpe_fit <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))

  list(
    fixef = lme4::fixef(fit),
    logLik = as.numeric(stats::logLik(fit)),
    sigma = stats::sigma(fit),
    vcov = vc$vcov,
    grp = vc$grp
  )
}

expect_mlpe_equivalent <- function(current, legacy, tolerance = 1e-7) {
  cur <- summarize_mlpe_fit(current)
  old <- summarize_mlpe_fit(legacy)

  expect_equal(cur$fixef, old$fixef, tolerance = tolerance)
  expect_equal(cur$logLik, old$logLik, tolerance = tolerance)
  expect_equal(cur$sigma, old$sigma, tolerance = tolerance)
  expect_equal(cur$vcov, old$vcov, tolerance = tolerance)
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

test_that("mlpe_rga accepts supplied ZZ together with keep", {
  n   <- 6
  df  <- make_mlpe_data(n)
  id  <- ResistanceGA2::To.From.ID(n)
  zz  <- ResistanceGA2:::ZZ.mat(id)
  keep <- rep(c(1L, 0L), length.out = nrow(df))

  expect_no_error(
    ResistanceGA2::mlpe_rga(
      y ~ x + (1 | pop),
      data = df,
      ZZ = zz,
      keep = keep
    )
  )
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

test_that("mlpe_data builds long-form dyad data from symmetric matrices", {
  y <- matrix(1:16, 4, 4)
  x <- matrix(seq(10, 25), 4, 4)

  dat <- ResistanceGA2::mlpe_data(response = y, x = x)

  expect_named(dat, c("from", "to", "x", "response"))
  expect_equal(nrow(dat), choose(4, 2))
  expect_equal(dat$response, ResistanceGA2::lower(y))
  expect_equal(dat$x, ResistanceGA2::lower(x))
})

test_that("mlpe_data supports explicit pair tables and keep filters", {
  pairs <- data.frame(from = c("a", "a", "b"), to = c("b", "c", "c"))

  dat <- ResistanceGA2::mlpe_data(
    response = c(1, 2, 3),
    x = c(4, 5, 6),
    pairs = pairs,
    keep = c(1L, 0L, 1L)
  )

  expect_equal(nrow(dat), 2L)
  expect_equal(dat$from, c("a", "b"))
  expect_equal(dat$to, c("b", "c"))
  expect_equal(dat$response, c(1, 3))
})

test_that("mlpe fits a dyadic random effect from endpoint columns", {
  y <- matrix(rnorm(25), 5)
  x <- matrix(rnorm(25), 5)
  dat <- ResistanceGA2::mlpe_data(response = y, x = x)

  fit <- ResistanceGA2::mlpe(
    response ~ x + (1 | pair),
    data = dat,
    pairs = c("from", "to")
  )

  expect_s4_class(fit, "lmerMod")
  expect_named(lme4::fixef(fit), c("(Intercept)", "x"))
})

test_that("mlpe preserves ordinary random intercepts alongside dyad terms", {
  set.seed(123)
  y <- matrix(rnorm(36), 6)
  x <- matrix(rnorm(36), 6)
  dat <- ResistanceGA2::mlpe_data(response = y, x = x)
  dat$species <- factor(rep(c("sp1", "sp2"), length.out = nrow(dat)))

  fit <- ResistanceGA2::mlpe(
    response ~ x + (1 | pair) + (1 | species),
    data = dat,
    pairs = c("from", "to")
  )

  expect_s4_class(fit, "lmerMod")
})

test_that("mlpe handles omitted pair observations via keep", {
  y <- matrix(rnorm(36), 6)
  x <- matrix(rnorm(36), 6)
  dat <- ResistanceGA2::mlpe_data(response = y, x = x)
  keep <- rep(c(TRUE, FALSE), length.out = nrow(dat))

  fit <- ResistanceGA2::mlpe(
    response ~ x + (1 | pair),
    data = dat,
    pairs = c("from", "to"),
    keep = keep
  )

  expect_s4_class(fit, "lmerMod")
})

test_that("mlpe and mlpe_rga match frozen legacy Gaussian fits", {
  set.seed(20260401)
  n <- 7
  y <- matrix(rnorm(n^2), n)
  x <- matrix(rnorm(n^2), n)
  id <- ResistanceGA2::To.From.ID(n)

  legacy_df <- data.frame(
    y = ResistanceGA2::lower(y),
    x = ResistanceGA2::lower(x),
    pop = id$pop1
  )
  long_df <- ResistanceGA2::mlpe_data(response = y, x = x)

  legacy_fit <- legacy_mlpe_rga(y ~ x + (1 | pop),
                                data = legacy_df,
                                REML = FALSE)
  current_rga <- ResistanceGA2::mlpe_rga(y ~ x + (1 | pop),
                                         data = legacy_df,
                                         REML = FALSE)
  current_mlpe <- ResistanceGA2::mlpe(response ~ x + (1 | pair),
                                      data = long_df,
                                      pairs = c("from", "to"),
                                      REML = FALSE)

  expect_mlpe_equivalent(current_rga, legacy_fit)
  expect_mlpe_equivalent(current_mlpe, legacy_fit)
})

test_that("mlpe and mlpe_rga match frozen legacy fits with omitted pairs", {
  set.seed(20260401)
  n <- 7
  y <- matrix(rnorm(n^2), n)
  x <- matrix(rnorm(n^2), n)
  id <- ResistanceGA2::To.From.ID(n)
  keep <- rep(c(1L, 0L), length.out = choose(n, 2))

  legacy_df <- data.frame(
    y = ResistanceGA2::lower(y),
    x = ResistanceGA2::lower(x),
    pop = id$pop1
  )
  long_df <- ResistanceGA2::mlpe_data(response = y, x = x)

  legacy_fit <- legacy_mlpe_rga(y ~ x + (1 | pop),
                                data = legacy_df,
                                REML = FALSE,
                                keep = keep)
  current_rga <- ResistanceGA2::mlpe_rga(y ~ x + (1 | pop),
                                         data = legacy_df,
                                         REML = FALSE,
                                         keep = keep)
  current_mlpe <- ResistanceGA2::mlpe(response ~ x + (1 | pair),
                                      data = long_df,
                                      pairs = c("from", "to"),
                                      REML = FALSE,
                                      keep = keep)

  expect_mlpe_equivalent(current_rga, legacy_fit)
  expect_mlpe_equivalent(current_mlpe, legacy_fit)
})

test_that("mlpe and mlpe_rga match frozen legacy Poisson fits", {
  set.seed(42)
  n <- 7
  y <- matrix(rpois(n^2, lambda = 4), n)
  x <- matrix(rnorm(n^2), n)
  id <- ResistanceGA2::To.From.ID(n)

  legacy_df <- data.frame(
    y = ResistanceGA2::lower(y),
    x = ResistanceGA2::lower(x),
    pop = id$pop1
  )
  long_df <- ResistanceGA2::mlpe_data(response = y, x = x)

  legacy_fit <- legacy_mlpe_rga(y ~ x + (1 | pop),
                                data = legacy_df,
                                family = poisson)
  current_rga <- ResistanceGA2::mlpe_rga(y ~ x + (1 | pop),
                                         data = legacy_df,
                                         family = poisson)
  current_mlpe <- ResistanceGA2::mlpe(response ~ x + (1 | pair),
                                      data = long_df,
                                      pairs = c("from", "to"),
                                      family = poisson)

  expect_mlpe_equivalent(current_rga, legacy_fit)
  expect_mlpe_equivalent(current_mlpe, legacy_fit)
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

test_that("MLPE.lmm matches frozen legacy outputs", {
  fixture <- make_mlpe_fixture()
  on.exit(unlink(fixture$file), add = TRUE)

  legacy_mat <- legacy_MLPE_lmm(fixture$matrix, fixture$genetic, REML = FALSE, scale = TRUE)
  current_mat <- ResistanceGA2::MLPE.lmm(fixture$matrix, fixture$genetic, REML = FALSE, scale = TRUE)
  expect_mlpe_equivalent(current_mat, legacy_mat)

  legacy_vec <- legacy_MLPE_lmm(fixture$vector, fixture$genetic, REML = FALSE, scale = FALSE)
  current_vec <- ResistanceGA2::MLPE.lmm(fixture$vector, fixture$genetic, REML = FALSE, scale = FALSE)
  expect_mlpe_equivalent(current_vec, legacy_vec)
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
