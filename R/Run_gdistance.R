#' Run gdistance cost/commute distance calculation
#'
#' Calculates pairwise cost or commute distances across a resistance surface
#' using the \code{gdistance} package.
#'
#' @param gdist.inputs Object created by \code{\link[ResistanceGA2]{gdist.prep}}.
#' @param r A \code{SpatRaster} object (from the \pkg{terra} package) or the
#'   path to a raster file readable by \code{terra::rast()}.
#' @param scl Logical. Scale the correction values (default \code{TRUE}). Set
#'   to \code{FALSE} to obtain absolute distance values. See
#'   \code{\link[gdistance]{geoCorrection}} for details.
#'
#' @return A numeric distance vector (or matrix) of pairwise cost/commute
#'   distances. Returns \code{-99999} on error or warning.
#'
#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r <- rast(nrows = 50, ncols = 50, vals = runif(2500))
#' coords <- matrix(runif(20, 0, 50), ncol = 2)
#' gdist.inputs <- gdist.prep(n.Pops = 10, samples = coords)
#' cd <- Run_gdistance(gdist.inputs, r)
#' }
Run_gdistance <- function(gdist.inputs,
                          r,
                          scl = TRUE) {
  out <- tryCatch(
    {
      if (!inherits(r, "SpatRaster")) {
        r <- terra::rast(r)
      }

      tr <- gdistance::transition(
        x = r,
        transitionFunction = gdist.inputs$transitionFunction,
        directions = gdist.inputs$directions
      )

      samples <- gdist.inputs$samples  # coordinate matrix

      if (is.null(gdist.inputs$keep)) {
        # All pairwise distances
        if (gdist.inputs$method == 'costDistance') {
          trC <- gdistance::geoCorrection(tr, "c", scl = scl)
          ret <- gdistance::costDistance(trC, samples)
          rm(trC)
        } else {
          trR <- gdistance::geoCorrection(tr, "r", scl = scl)
          ret <- gdistance::commuteDistance(trR, samples) / 1000
          rm(trR)
        }

      } else {
        # Run on selected pairs only
        kp  <- which(gdist.inputs$keep == 1)
        ret <- vector(mode = "numeric", length = length(kp))
        id  <- as.matrix(gdist.inputs$ID)
        id  <- apply(id, 2, as.numeric)

        if (gdist.inputs$method == 'costDistance') {
          trC <- gdistance::geoCorrection(tr, "c", scl = scl)
          for (i in seq_along(kp)) {
            pts      <- as.integer(id[kp[i], ])
            ret[[i]] <- gdistance::costDistance(trC, samples[pts, , drop = FALSE])
          }
          rm(trC)
        } else {
          trR <- gdistance::geoCorrection(tr, "r", scl = scl)
          for (i in seq_along(kp)) {
            pts      <- as.integer(id[kp[i], ])
            ret[[i]] <- gdistance::commuteDistance(trR, samples[pts, , drop = FALSE]) / 1000
          }
          rm(trR)
        }
      }

      rm(tr, r)
      gc()
      return(ret)
    },
    warning = function(w) {
      return(-99999)
    },
    error = function(e) {
      return(-99999)
    }
  )
  return(out)
}
