% Generated manually for site-level pairwise environmental helpers
\name{site_pairwise_dist}
\alias{site_pairwise_dist}
\title{Build Pairwise Environmental-Distance Covariates}
\usage{
site_pairwise_dist(data, vars = NULL, scale = TRUE, method = "euclidean",
  keep = NULL, pop_n = NULL, name = "env_dist", ...)
}
\arguments{
\item{data}{Numeric vector, matrix, or data frame of site-level environmental variables.}

\item{vars}{Optional character or integer index selecting columns from \code{data}.}

\item{scale}{Logical. If \code{TRUE}, standardize each site-level variable before calculating environmental distance. Constant variables are converted to zeros.}

\item{method}{Distance metric passed to \code{\link[stats]{dist}}.}

\item{keep}{Optional logical/integer vector used to retain a subset of the final pairwise observations.}

\item{pop_n}{Optional integer vector giving the number of sampled individuals per population. When supplied, population-level environmental distances are expanded to the individual-level pair structure used by \code{\link[ResistanceGA2]{jl.prep}} with \code{pop2ind}.}

\item{name}{Name of the output covariate column. Default = \code{"env_dist"}.}

\item{...}{Additional arguments passed to \code{\link[stats]{dist}}, such as \code{p} when \code{method = "minkowski"}.}
}
\value{
A one-column data frame containing pairwise environmental distance.
}
\description{
Convert one or more site-level environmental variables into a single pairwise environmental-distance covariate aligned with the lower-triangle ordering used throughout \pkg{ResistanceGA2}.
}
\examples{
env <- data.frame(
  temp = c(10, 12, 15, 11),
  precip = c(100, 120, 80, 110)
)

site_pairwise_dist(env)
site_pairwise_dist(env, method = "manhattan", name = "climate_dist")
}
