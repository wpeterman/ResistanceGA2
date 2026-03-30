#' Optimize a single resistance surface (internal fitness function)
#'
#' Transforms and evaluates one resistance surface during GA optimization.
#' Designed to be called from \code{SS_optim}.
#'
#' @param PARM Parameter vector. For categorical surfaces: one value per class
#'   (resistance multipliers). For continuous surfaces: transformation code,
#'   shape, and maximum value.
#' @param Resistance \code{SpatRaster} resistance surface to optimize.
#' @param gdist.inputs Object from \code{\link{gdist.prep}}.
#' @param jl.inputs Object from \code{\link{jl.prep}}.
#' @param GA.inputs Object from \code{\link{GA.prep}}.
#' @param Min.Max \code{'min'} or \code{'max'} (default). Direction of
#'   optimization passed by the GA.
#' @param iter Integer counter indicating which surface (layer) is being
#'   optimized.
#' @param quiet Logical. If \code{TRUE}, suppress per-iteration console output.
#'   Default = \code{FALSE}.
#'
#' @return Scalar objective function value (negative AIC, R², or log-likelihood).
#' @noRd
#' @author Bill Peterman <Peterman.73@@osu.edu>
Resistance.Opt_single <-
  function(PARM,
           Resistance,
           gdist.inputs = NULL,
           jl.inputs    = NULL,
           GA.inputs,
           Min.Max      = "max",
           iter         = NULL,
           quiet        = FALSE) {

    t1   <- Sys.time()
    iter <- iter %||% 1L

    method       <- GA.inputs$method
    select.trans <- GA.inputs$select.trans
    r            <- Resistance
    keep         <- 1L

    # Categorical surface -------------------------------------------------------
    if (GA.inputs$surface.type[iter] == "cat") {
      PARM <- PARM / min(PARM)
      if (max(PARM) > GA.inputs$max.cat) {
        PARM <- SCALE.vector(PARM, 1, GA.inputs$max.cat)
      }
      lev <- terra::unique(r)[[1]]
      df  <- data.frame(id = lev, PARM = PARM)
      r   <- terra::subst(r, from = df$id, to = df$PARM)
      equation <- NA_integer_   # not applicable for categorical

    } else {
      # Continuous surface ------------------------------------------------------
      r        <- SCALE(r, 0, 10)
      equation <- floor(PARM[1])
      SHAPE    <- PARM[2]

      if (equation %in% select.trans[[iter]]) {
        rick.eq <- equation %in% c(2L, 4L, 6L, 8L)
        if (rick.eq && SHAPE > 5) {
          equation <- 9L
          keep     <- 0L
        }

        r <- switch(
          as.character(equation),
          "1" = Inv.Rev.Monomolecular(r, PARM),
          "2" = Inv.Rev.Ricker(r, PARM),
          "3" = Monomolecular(r, PARM),
          "4" = Ricker(r, PARM),
          "5" = Rev.Monomolecular(r, PARM),
          "6" = Rev.Ricker(r, PARM),
          "7" = Inv.Monomolecular(r, PARM),
          "8" = Inv.Ricker(r, PARM),
          { keep <- 0L; (r * 0) + 1 }   # Distance fallback
        )
      }
    }

    # Evaluate objective function -----------------------------------------------
    surface_ok <- (GA.inputs$surface.type[iter] == "cat") ||
                  (!is.na(equation) && equation %in% select.trans[[iter]])

    if (!surface_ok || keep == 0L) {
      obj.func.opt <- -99999
    } else {
      rclmat <- matrix(c(1e-100, 1e-6, 1e-6,
                         1e-6,   1e6,  NA,        # keep values in range (pass-through handled by classify)
                         1e6,    Inf,  1e6,
                         -1,     0,    1e6),
                       ncol = 3, byrow = TRUE)

      # Reclassify: values <= 0 → 1e6; values > 1e6 → 1e6; tiny values → 1e-6
      rclmat2 <- matrix(c(-Inf, 1e-100, 1e-6,
                          1e6,  Inf,    1e6),
                        ncol = 3, byrow = TRUE)
      r <- terra::classify(r, rclmat2, include.lowest = TRUE)

      obj.func.opt <- -99999   # default (failure)

      # gdistance ---------------------------------------------------------------
      if (!is.null(gdist.inputs)) {
        if (mean(terra::values(r), na.rm = TRUE) == 0) {
          obj.func.opt <- -99999
        } else {
          cd <- try(Run_gdistance(gdist.inputs, r), silent = TRUE)

          if (inherits(cd, "try-error") || identical(cd, -99999)) {
            obj.func.opt <- -99999
          } else {
            dat <- data.frame(gd  = gdist.inputs$response,
                              cd  = scale(c(cd)),
                              pop = gdist.inputs$ID$pop1)

            fit.mod <- mlpe_rga(formula = gd ~ cd + (1 | pop),
                                data    = dat,
                                ZZ      = gdist.inputs$ZZ,
                                REML    = FALSE)

            if (lme4::fixef(fit.mod)["cd"] < 0) {
              obj.func.opt <- -99999
            } else {
              obj.func.opt <- .obj_value(fit.mod, method)
            }
          }
        }
        rm(cd, r)
        gc()
      }

      # Julia / Circuitscape ----------------------------------------------------
      if (!is.null(jl.inputs)) {
        if (mean(terra::values(r), na.rm = TRUE) == 0) {
          obj.func.opt <- -99999
        } else {
          cd <- try(Run_CS.jl(jl.inputs, r), silent = TRUE)

          if (inherits(cd, "try-error") || identical(cd, -99999)) {
            obj.func.opt <- -99999
          } else {
            dat     <- jl.inputs$df
            dat$cd  <- scale(cd)

            fit.mod <- mlpe_rga(formula = gd ~ cd + (1 | pop),
                                data    = dat,
                                ZZ      = jl.inputs$ZZ,
                                REML    = FALSE)

            if (lme4::fixef(fit.mod)["cd"] < 0) {
              obj.func.opt <- -99999
            } else {
              obj.func.opt <- .obj_value(fit.mod, method)
            }
          }
        }
      }
    }

    rt <- Sys.time() - t1
    if (!quiet) {
      cat(paste0("\t", "Iteration took ", round(rt, 2), " seconds\n"))
      cat(paste0("\t", method, " = ", round(obj.func.opt, 4), "\n\n"))
      if (!is.null(iter) && GA.inputs$surface.type[iter] != "cat" && !is.na(equation)) {
        EQ <- .eq_name(equation)
        cat(paste0("\t", EQ, " | Shape = ", PARM[2], " | Max = ", PARM[3], "\n\n"))
      }
    }

    gc()
    if (!is.null(GA.inputs$opt.digits)) {
      return(round(obj.func.opt, GA.inputs$opt.digits))
    }
    return(obj.func.opt)
  }

# Internal helpers ------------------------------------------------------------

.obj_value <- function(fit.mod, method) {
  switch(method,
    "AIC" = suppressWarnings(AIC(fit.mod)) * -1,
    "R2"  = suppressWarnings(MuMIn::r.squaredGLMM(fit.mod))[[1]],
    suppressWarnings(logLik(fit.mod))[[1]]
  )
}

.eq_name <- function(equation) {
  switch(as.character(equation),
    "1" = "Inverse-Reverse Monomolecular",
    "2" = "Inverse-Reverse Ricker",
    "3" = "Monomolecular",
    "4" = "Ricker",
    "5" = "Reverse Monomolecular",
    "6" = "Reverse Ricker",
    "7" = "Inverse Monomolecular",
    "8" = "Inverse Ricker",
    "Distance"
  )
}
