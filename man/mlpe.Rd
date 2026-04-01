% Generated manually for new MLPE interface
\name{mlpe}
\alias{mlpe}
\title{Fit MLPE Models With an lme4-Style Formula}
\usage{
mlpe(formula, data, pairs, REML = FALSE, keep = NULL, ...)
}
\arguments{
\item{formula}{Mixed-effects model formula.}

\item{data}{Data frame containing the response, covariates, and endpoint columns.}

\item{pairs}{Definition of one or more MLPE dyad terms. Supply either a character vector of length two, which is interpreted as \code{list(pair = c("from", "to"))}, or a named list mapping random-effect term names in \code{formula} to the two endpoint columns used to construct each dyadic random intercept.}

\item{REML}{Logical. If TRUE, fit by restricted maximum likelihood for linear mixed models.}

\item{keep}{Optional logical/integer vector used to retain a subset of rows from \code{data} before fitting.}

\item{...}{Additional arguments passed to \code{\link[lme4]{glFormula}} when fitting generalized models, such as \code{family}.}
}
\value{
A fitted \code{merMod} object.
}
\description{
Fit maximum-likelihood population effects (MLPE) models using standard \code{lme4}-style formulas while defining dyadic random-intercept terms from explicit endpoint columns.
}
\details{
Dyadic random effects are specified in \code{formula} using standard random intercept syntax, for example \code{(1 | pair)}. The corresponding endpoint columns are declared through \code{pairs}, e.g. \code{pairs = c("from", "to")} or \code{pairs = list(pair = c("from", "to"))}.

Additional ordinary \code{lme4} random effects may be included in the same model. The dyadic covariance contribution is rebuilt after any row subsetting so omitted pairwise observations are handled correctly.
}
\examples{
y <- matrix(rnorm(25), 5)
x <- matrix(rnorm(25), 5)

dat <- mlpe_data(response = y, resistance = x)

fit <- mlpe(
  response ~ resistance + (1 | pair),
  data = dat,
  pairs = c("from", "to")
)
}
