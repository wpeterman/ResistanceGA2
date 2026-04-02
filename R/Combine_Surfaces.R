#' Combine multiple resistance surfaces into a composite surface
#'
#' Combines multiple resistance surfaces by transforming each layer with
#' specified parameters and summing the results.
#'
#' @param PARM Parameters to transform continuous or categorical resistance
#'   surfaces. Provide a vector with parameters in the order of resistance
#'   surfaces as output by \code{\link{MS_optim}}.
#' @param gdist.inputs Object from \code{\link{gdist.prep}}. Supply when
#'   optimizing with gdistance.
#' @param jl.inputs Object from \code{\link{jl.prep}}. Supply when optimizing
#'   with Circuitscape via Julia.
#' @param GA.inputs Object from \code{\link{GA.prep}}.
#' @param out Directory to write the combined \code{.asc} file. Default =
#'   \code{NULL} (no file written).
#' @param File.name Name for the output \code{.asc} file. Default is the
#'   surface names joined by \code{"."}.
#' @param rescale Logical. If \code{TRUE} (default), the combined surface is
#'   divided by its minimum value so the minimum resistance equals 1.
#' @param p.contribution Logical. If \code{TRUE}, returns a list with (1) the
#'   combined raster surface and (2) the mean proportional contribution of each
#'   input surface to the combined resistance values.
#'
#' @return A \code{SpatRaster} of the combined resistance surface (sum of all
#'   transformed/reclassified layers), or a list when
#'   \code{p.contribution = TRUE}.
#'
#' @details
#' For continuous surfaces, \code{PARM} should specify, per layer, three
#' values: (1) transformation code (1–9), (2) shape, and (3) maximum value.
#' Transformation codes:
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
#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @examples
#' \dontrun{
#' # After running MS_optim(), use the best-fit PARM vector:
#' combined <- Combine_Surfaces(
#'   PARM        = ms.out$PARM,
#'   gdist.inputs = gdist.inputs,
#'   GA.inputs   = GA.inputs
#' )
#' }

Combine_Surfaces <-
  function(PARM,
           gdist.inputs = NULL,
           jl.inputs    = NULL,
           GA.inputs,
           out          = NULL,
           File.name    = paste(GA.inputs$parm.type$name, collapse = "."),
           rescale      = TRUE,
           p.contribution = FALSE) {

    if (!is.null(gdist.inputs)) {
      EXPORT.dir <- out
    }

    if (!is.null(jl.inputs)) {
      EXPORT.dir <- out
    }

    materialize_raster <- function(x) {
      if (inherits(x, "PackedSpatRaster")) {
        terra::unwrap(x)
      } else {
        x
      }
    }

    select.trans <- GA.inputs$select.trans
    r    <- materialize_raster(GA.inputs$Resistance.stack)
    keep <- 1

    for (i in seq_len(GA.inputs$n.layers)) {
      if (GA.inputs$surface.type[i] == "cat") {
        parm <- PARM[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]

        if (is.na(sum(parm))) {
          parm <- rep(1, length(parm))
          keep <- 0
        }

        parm <- parm / min(parm)
        if (max(parm) > GA.inputs$max.cat) {
          parm <- SCALE.vector(parm, 1, GA.inputs$max.cat)
        }

        lev <- terra::unique(r[[i]])[[1]]
        df  <- data.frame(id = lev, parm = parm)
        r[[i]] <- terra::subst(r[[i]], from = df$id, to = df$parm)

      } else {
        rast <- SCALE(data = r[[i]], MIN = 0, MAX = 10)
        parm <- PARM[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]

        if (is.na(parm[1])) { parm[1] <- 9; keep <- 0 }
        if (is.na(parm[2])) { parm[2] <- 1; parm[1] <- 9; keep <- 0 }
        if (is.na(parm[3])) { parm[3] <- 2; parm[1] <- 9; keep <- 0 }

        equation   <- floor(parm[1])
        SHAPE      <- parm[2]
        rick.eq    <- equation %in% c(2, 4, 6, 8)

        if (rick.eq && SHAPE > 5) {
          equation <- 9
          keep     <- 0
        }

        if (equation %in% select.trans[[i]] && keep == 1) {
          keep <- 1
        } else {
          equation <- 9
          keep     <- 0
        }

        r[[i]] <- switch(
          as.character(equation),
          "1" = Inv.Rev.Monomolecular(rast, parm),
          "2" = Inv.Rev.Ricker(rast, parm),
          "3" = Monomolecular(rast, parm),
          "4" = Ricker(rast, parm),
          "5" = Rev.Monomolecular(rast, parm),
          "6" = Rev.Ricker(rast, parm),
          "7" = Inv.Monomolecular(rast, parm),
          "8" = Inv.Ricker(rast, parm),
          (rast * 0) + 1  # Distance / unused
        )
      }
    }

    multi_surface <- sum(r)

    if (keep == 0) {
      multi_surface <- r[[1]] * 0
    }

    # Helper: check for all-NA surface
    all_na <- function(surf) {
      sum(!is.na(terra::values(surf))) == 0
    }

    # Rescale / cap helpers
    rescale_surface <- function(surf) {
      mn <- terra::global(surf, "min", na.rm = TRUE)[[1]]
      surf / mn
    }
    cap_surface <- function(surf) {
      if (terra::global(surf, "max", na.rm = TRUE)[[1]] > 1e6) {
        SCALE(surf, 1, 1e6)
      } else {
        surf
      }
    }

    if (p.contribution) {
      cont.list <- vector("list", GA.inputs$n.layers)
      for (i in seq_len(GA.inputs$n.layers)) {
        p.cont       <- r[[i]] / multi_surface
        mean.cont    <- terra::global(p.cont, "mean", na.rm = TRUE)[[1]]
        cont.list[[i]] <- data.frame(surface = GA.inputs$layer.names[i],
                                     mean    = mean.cont)
      }

      if (keep != 0 && rescale)   multi_surface <- rescale_surface(multi_surface)
      if (keep != 0)              multi_surface <- cap_surface(multi_surface)

      if (all_na(multi_surface)) {
        keep          <- 0
        multi_surface <- r[[1]] * 0
      }

      cont.df <- plyr::ldply(cont.list)
      return(list(percent.contribution = cont.df,
                  combined.surface     = multi_surface))
    }

    # p.contribution == FALSE
    if (keep != 0 && rescale) multi_surface <- rescale_surface(multi_surface)
    if (keep != 0)            multi_surface <- cap_surface(multi_surface)

    if (all_na(multi_surface)) {
      keep          <- 0
      multi_surface <- r[[1]] * 0
    }

    if (!is.null(out)) {
      terra::writeRaster(
        x        = multi_surface,
        filename = file.path(out, paste0(File.name, ".asc")),
        overwrite = TRUE
      )
    }

    rm(r)
    gc()
    return(multi_surface)
  }
