#' Optimize multiple resistance surfaces with covariate formula (internal fitness function)
#'
#' Creates and evaluates a composite resistance surface during GA optimization,
#' using user-supplied formulas for the MLPE model (supports covariates).
#' Designed to be called from \code{MS_optim}.
#'
#' @param PARM Parameter vector for all surfaces in the order of layers in
#'   \code{GA.inputs$Resistance.stack}.
#' @param gdist.inputs Object from \code{\link{gdist.prep}}.
#' @param jl.inputs Object from \code{\link{jl.prep}}.
#' @param GA.inputs Object from \code{\link{GA.prep}}.
#' @param Min.Max \code{'min'} or \code{'max'} (default).
#' @param quiet Logical. Suppress per-iteration console output? Default =
#'   \code{FALSE}.
#'
#' @return Scalar objective function value.
#' @noRd
#' @author Bill Peterman <Peterman.73@@osu.edu>
Resistance.Opt_multi.cov <- function(PARM,
                                     gdist.inputs = NULL,
                                     jl.inputs    = NULL,
                                     GA.inputs,
                                     Min.Max      = "max",
                                     quiet        = FALSE) {
  materialize_raster <- function(x) {
    if (inherits(x, "PackedSpatRaster")) {
      terra::unwrap(x)
    } else {
      x
    }
  }

  t1        <- proc.time()[3]
  method    <- GA.inputs$method
  File.name <- "resist_surface"
  worker.inputs <- GA.inputs
  worker.inputs$Resistance.stack <- materialize_raster(GA.inputs$Resistance.stack)

  materialize_raster <- function(x) {
    if (inherits(x, "PackedSpatRaster")) {
      terra::unwrap(x)
    } else {
      x
    }
  }

  worker.inputs <- GA.inputs
  worker.inputs$Resistance.stack <-
    materialize_raster(GA.inputs$Resistance.stack)

  obj.func.opt <- -99999

  # gdistance -----------------------------------------------------------------
  if (!is.null(gdist.inputs)) {
    r <- Combine_Surfaces(
      PARM         = PARM,
      gdist.inputs = gdist.inputs,
      GA.inputs    = worker.inputs,
      out          = NULL,
      File.name    = File.name,
      rescale      = FALSE
    )

    if (mean(terra::values(r), na.rm = TRUE) != 0) {
      cd <- try(
        Run_gdistance(gdist.inputs, r, return.error.value = TRUE),
        silent = TRUE
      )

      if (!inherits(cd, "try-error") && !identical(cd, -99999)) {
        dat    <- gdist.inputs$df
        dat$cd <- scale(c(cd))

        fit.mod <- mlpe_rga(formula = gdist.inputs$formula,
                            data    = dat,
                            ZZ      = gdist.inputs$ZZ,
                            REML    = FALSE)

        if (lme4::fixef(fit.mod)["cd"] >= 0) {
          obj.func.opt <- .obj_value(fit.mod, method)
        }
      }
    }
    rm(r)
    gc()
  }

  # Julia / Circuitscape ------------------------------------------------------
  if (!is.null(jl.inputs)) {
    r <- Combine_Surfaces(
      PARM      = PARM,
      jl.inputs = jl.inputs,
      GA.inputs = worker.inputs,
      out       = NULL,
      File.name = File.name,
      rescale   = FALSE
    )

    if (mean(terra::values(r), na.rm = TRUE) != 0) {
      cd <- try(Run_CS.jl(jl.inputs, r, full.mat = FALSE), silent = TRUE)

      if (!inherits(cd, "try-error") && !identical(cd, -99999)) {
        dat <- jl.inputs$df

        if (is.null(nrow(cd))) {
          dat$cd <- scale(cd)
        } else {
          dat$cd <- scale(lower(cd))
        }

        fit.mod <- mlpe_rga(formula = jl.inputs$formula,
                            data    = dat,
                            ZZ      = jl.inputs$ZZ,
                            REML    = FALSE)

        if (lme4::fixef(fit.mod)["cd"] >= 0) {
          obj.func.opt <- .obj_value(fit.mod, method)
        }
      }
    }
  }

  rt <- proc.time()[3] - t1
  if (!quiet) {
    cat(paste0("\t", "Iteration took ", round(rt, 2), " seconds\n"))
    cat(paste0("\t", method, " = ", round(obj.func.opt, 4), "\n\n"))
  }
  gc()

  if (!is.null(GA.inputs$opt.digits)) {
    return(round(obj.func.opt, GA.inputs$opt.digits))
  }
  return(obj.func.opt)
}
