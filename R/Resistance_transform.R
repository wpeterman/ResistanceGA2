#' Apply a transformation to a continuous resistance surface
#'
#' Applies one of eight resistance transformations (or a distance/flat
#' transformation) to a continuous resistance surface.
#'
#' @param transformation Transformation to apply. Either a numeric code (1–9)
#'   or the name of the transformation (see Details).
#' @param shape Shape parameter value.
#' @param max Maximum value parameter.
#' @param r Resistance surface as a single-layer \code{SpatRaster} object.
#' @param out Directory path for exporting the transformed surface as an
#'   \code{.asc} file. Default = \code{NULL} (no file written).
#'
#' @return A \code{SpatRaster} with the transformation applied.
#'
#' @details
#' Valid values for \code{transformation}:
#' \tabular{ll}{
#'   \tab 1 = "Inverse-Reverse Monomolecular"\cr
#'   \tab 2 = "Inverse-Reverse Ricker"\cr
#'   \tab 3 = "Monomolecular"\cr
#'   \tab 4 = "Ricker"\cr
#'   \tab 5 = "Reverse Monomolecular"\cr
#'   \tab 6 = "Reverse Ricker"\cr
#'   \tab 7 = "Inverse Monomolecular"\cr
#'   \tab 8 = "Inverse Ricker"\cr
#'   \tab 9 = "Distance"\cr
#' }
#'
#' The Distance transformation sets all cell values equal to 1. Because the
#' Ricker function can approximate a monomolecular shape at high shape
#' parameter values, whenever a shape parameter > 6 is selected with a Ricker
#' family transformation, the transformation automatically reverts to Distance.
#' In most optimization workflows, monomolecular families should be treated as
#' the recommended default, and Ricker families should only be explored when
#' there is strong biological or ecological support for a unimodal or
#' quadratic-type relationship. Otherwise, the additional flexibility may
#' increase the risk of overfitting.
#'
#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @examples
#' r_trans <- Resistance.tran(
#'   transformation = "Monomolecular",
#'   shape = 2.5,
#'   max = 100,
#'   r = raster_orig[["cont_orig"]]
#' )
#'
#' r_trans

Resistance.tran <- function(transformation,
                            shape,
                            max,
                            r,
                            out = NULL) {

  R <- .validate_spatraster(r, arg = "r", nlyr = 1L)
  NAME <- names(R)

  if (is.numeric(transformation)) {
    parm <- c(transformation, shape, max)
  } else {
    parm <- c(get.EQ(transformation), shape, max)
  }

  R <- SCALE(data = R, MIN = 0, MAX = 10)

  equation <- floor(parm[1])  # parameter can range 1–9.99

  r <- switch(
    as.character(equation),
    "1" = { Inv.Rev.Monomolecular(R, parm) },
    "2" = { Inv.Rev.Ricker(R, parm)        },
    "3" = { Monomolecular(R, parm)         },
    "4" = { Ricker(R, parm)                },
    "5" = { Rev.Monomolecular(R, parm)     },
    "6" = { Rev.Ricker(R, parm)            },
    "7" = { Inv.Monomolecular(R, parm)     },
    "8" = { Inv.Ricker(R, parm)            },
    (R * 0) + 1  # Distance (equation == 9 or unknown)
  )

  if (terra::global(r, "max", na.rm = TRUE)[[1]] > 1e6) {
    r <- SCALE(r, 1, 1e6)
  }

  if (!is.null(out)) {
    terra::writeRaster(
      x        = r,
      filename = file.path(out, paste0(NAME, ".asc")),
      overwrite = TRUE
    )
  }

  return(r)
}
