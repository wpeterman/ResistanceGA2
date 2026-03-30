#' Prepare data for optimization using \code{gdistance}
#'
#' Creates a necessary input for optimizing resistance surfaces based on
#' pairwise cost distances, implemented using the \code{gdistance} library.
#'
#' @param n.Pops The number of populations that are being assessed.
#' @param response Vector of pairwise genetic distances (lower half of pairwise
#'   matrix).
#' @param samples Either the path to a tab-delimited .txt file containing xy
#'   coordinates (columns 2 and 3, first column ignored), a two-column matrix
#'   with x values in column 1 and y values in column 2, or a
#'   \code{\link[terra]{SpatVector}} of points.
#' @param covariates Data frame of additional covariates to include in the MLPE
#'   model during optimization.
#' @param formula If covariates are included in the model, specify the R formula
#'   for the fixed effects portion of the MLPE model (e.g.,
#'   \code{response ~ covariate}). The \code{response} term will use the values
#'   supplied to the \code{response} parameter; covariate names must match
#'   column names in \code{covariates}.
#' @param transitionFunction Function to calculate the \code{gdistance}
#'   TransitionLayer. See \code{\link[gdistance]{transition}}.
#'   Default = \code{function(x) 1 / mean(x)}.
#' @param directions Directions in which cells are connected (4, 8, 16, or
#'   other). Default = 8.
#' @param longlat Logical. If \code{TRUE}, a \code{\link[gdistance]{geoCorrection}}
#'   will be applied to the transition matrix. Default = \code{FALSE}.
#' @param method Pairwise distance method: \code{'commuteDistance'} (default,
#'   equivalent to Circuitscape resistance distance) or \code{'costDistance'}
#'   (least-cost path distance). See \code{\link[gdistance]{costDistance}} and
#'   \code{\link[gdistance]{commuteDistance}}.
#' @param min.max_dist NOT YET SUPPORTED. Optional two-element vector
#'   \code{c(min, max)} specifying the Euclidean distance range for pairwise
#'   comparisons. Pairs outside this range are omitted.
#' @param keep NOT YET SUPPORTED. An optional vector of length equal to the
#'   number of pairwise observations, with 1 to keep an observation and 0 to
#'   drop it. Can be used in conjunction with, or in place of,
#'   \code{min.max_dist}.
#'
#' @return A named list of inputs required by the optimization functions.
#'
#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' set.seed(42)
#' coords <- matrix(runif(20, 0, 100), ncol = 2)
#' gd <- runif(choose(10, 2))
#'
#' gdist.inputs <- gdist.prep(
#'   n.Pops = 10,
#'   response = gd,
#'   samples = coords,
#'   method = "commuteDistance"
#' )
#' }

gdist.prep <-
  function(n.Pops,
           response = NULL,
           samples,
           covariates = NULL,
           formula = NULL,
           transitionFunction = function(x) 1 / mean(x),
           directions = 8,
           longlat = FALSE,
           method = 'commuteDistance',
           min.max_dist = NULL,
           keep = NULL) {

    if (method != 'commuteDistance') {
      method <- 'costDistance'
    }

    if (!is.null(response)) {
      if (!is.vector(response)) {
        stop("'response' must be a single-column vector of pairwise genetic distances.")
      }
    }

    # Validate covariates
    if (!is.null(covariates) && !is.data.frame(covariates)) {
      stop("'covariates' must be a data frame.")
    }
    if (!is.null(covariates) && !is.null(response) &&
        nrow(covariates) != length(response)) {
      stop("'response' and 'covariates' must have the same number of observations.")
    }

    # Parse sample locations into a coordinate matrix -------------------------
    if (is.matrix(samples)) {
      if (ncol(samples) > 2) {
        stop("The coordinate matrix has too many columns; supply x in column 1 and y in column 2.")
      }
      sp <- samples
    } else if (inherits(samples, "SpatVector")) {
      sp <- terra::crds(samples)
    } else if (is.data.frame(samples)) {
      sp <- as.matrix(samples)
    } else if (is.character(samples)) {
      if (!file.exists(samples)) {
        stop("The path to the samples file does not exist: ", samples)
      }
      sp <- as.matrix(read.delim(samples, header = FALSE)[, -1])
    } else {
      stop("'samples' must be a coordinate matrix, data frame, SpatVector, or file path.")
    }

    if (n.Pops != nrow(sp)) {
      stop("'n.Pops' (", n.Pops, ") does not equal the number of sample locations (",
           nrow(sp), ").")
    }

    # Build ID / ZZ matrices --------------------------------------------------
    if (!is.null(keep)) {
      ID <- To.From.ID(n.Pops)
      if (length(keep) != nrow(ID)) {
        stop("'keep' vector length must equal the number of pairwise combinations (",
             nrow(ID), ").")
      }
      ZZ <- ZZ.mat(ID, keep)
    } else {
      ID <- To.From.ID(n.Pops)
      ZZ <- ZZ.mat(ID)
    }

    # Build data frame and formula --------------------------------------------
    df <- NULL
    if (!is.null(response)) {
      if (!is.null(covariates)) {
        df <- data.frame(gd = response, covariates, pop = ID$pop1)
      } else {
        df <- data.frame(gd = response, pop = ID$pop1)
      }

      if (!is.null(formula)) {
        formula <- update(formula, gd ~ . + cd + (1 | pop))
      } else {
        formula <- gd ~ cd + (1 | pop)
      }

      # Warn if genetic distance decreases with Euclidean distance
      m <- lm(gd ~ c(dist(sp)), data = df)
      if (coef(m)[2] < 0) {
        warning(
          "Genetic distance decreases with Euclidean distance. ",
          "This is likely to result in a failed optimization.\n",
          "Check your measure carefully and consider subtracting values from 1 ",
          "to reverse the relationship."
        )
      }
    }

    ret <- list(
      response           = response,
      samples            = sp,
      covariates         = covariates,
      formula            = formula,
      transitionFunction = transitionFunction,
      directions         = directions,
      ID                 = ID,
      ZZ                 = ZZ,
      keep               = keep,
      n.Pops             = n.Pops,
      longlat            = longlat,
      method             = method,
      df                 = df
    )

    # Min-Max Distance (not yet supported) ------------------------------------
    if (!is.null(min.max_dist)) {
      stop("The 'min.max_dist' feature is not yet supported.")
    }

    return(ret)
  }
