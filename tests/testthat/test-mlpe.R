# Helper: generate a small MLPE-ready data frame
make_mlpe_data <- function(n = 5, seed = 1) {
  set.seed(seed)
  m_y <- matrix(rnorm(n^2), n, n)
  m_x <- matrix(rnorm(n^2), n, n)
  id  <- To.From.ID(n)
  data.frame(
    y   = lower(m_y),
    x   = lower(m_x),
    pop = id$pop1
  )
}

test_that("mlpe_rga returns a lmerMod object", {
  df  <- make_mlpe_data(5)
  fit <- mlpe_rga(y ~ x + (1 | pop), data = df)
  expect_s4_class(fit, "lmerMod")
})

test_that("mlpe_rga AIC is finite when REML = FALSE", {
  df  <- make_mlpe_data(5)
  fit <- mlpe_rga(y ~ x + (1 | pop), data = df, REML = FALSE)
  expect_true(is.finite(AIC(fit)))
})

test_that("mlpe_rga fixef has expected names", {
  df  <- make_mlpe_data(6)
  fit <- mlpe_rga(y ~ x + (1 | pop), data = df)
  fe  <- lme4::fixef(fit)
  expect_named(fe, c("(Intercept)", "x"))
})

test_that("mlpe_rga REML and ML give different log-likelihoods", {
  df       <- make_mlpe_data(6)
  fit_ml   <- mlpe_rga(y ~ x + (1 | pop), data = df, REML = FALSE)
  fit_reml <- mlpe_rga(y ~ x + (1 | pop), data = df, REML = TRUE)
  expect_false(isTRUE(all.equal(
    as.numeric(logLik(fit_ml)),
    as.numeric(logLik(fit_reml))
  )))
})

test_that("mlpe_rga accepts formula as character string", {
  df <- make_mlpe_data(5)
  expect_no_error(mlpe_rga("y ~ x + (1 | pop)", data = df))
})

test_that("mlpe_rga keep argument subsets observations", {
  df   <- make_mlpe_data(6)
  keep <- rep(c(1L, 0L), length.out = nrow(df))
  fit  <- mlpe_rga(y ~ x + (1 | pop), data = df, keep = keep)
  expect_s4_class(fit, "lmerMod")
})

test_that("mlpe_rga with supplied ZZ matches auto-computed ZZ", {
  n   <- 5
  df  <- make_mlpe_data(n)
  id  <- To.From.ID(n)
  ZZ  <- ResistanceGA2:::ZZ.mat(id)
  fit_auto <- mlpe_rga(y ~ x + (1 | pop), data = df)
  fit_man  <- mlpe_rga(y ~ x + (1 | pop), data = df, ZZ = ZZ)
  expect_equal(
    as.numeric(logLik(fit_auto)),
    as.numeric(logLik(fit_man))
  )
})

test_that("mlpe_rga Poisson glmer runs without error", {
  set.seed(42)
  n   <- 5
  id  <- To.From.ID(n)
  df  <- data.frame(
    y   = rpois(choose(n, 2), lambda = 5),
    x   = rnorm(choose(n, 2)),
    pop = id$pop1
  )
  expect_no_error(
    mlpe_rga(y ~ x + (1 | pop), data = df, family = poisson)
  )
})
