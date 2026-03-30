############ RUN JULIA ############

#' Run Circuitscape via Julia
#'
#' Executes Circuitscape through Julia (via \pkg{JuliaConnectoR}) and returns
#' pairwise resistance distances.
#'
#' @param jl.inputs Object from \code{\link{jl.prep}}. When provided, all other
#'   Julia/CS settings are taken from this object.
#' @param r A \code{SpatRaster} (from \pkg{terra}) or path to a raster file
#'   readable by \code{terra::rast()}.
#' @param CurrentMap Logical. Generate the cumulative current map?
#'   Default = \code{FALSE}.
#' @param full.mat Logical. Return the full square distance matrix instead of
#'   only the lower triangle? Default = \code{FALSE}. Cannot be requested when
#'   \code{pairs_to_include} is used.
#' @param EXPORT.dir Directory where Circuitscape results are written. Required
#'   when \code{CurrentMap = TRUE}.
#' @param output \code{"matrix"} (default) returns the lower-half pairwise
#'   resistance vector; \code{"raster"} returns a \code{SpatRaster} of the
#'   cumulative current map (requires \code{CurrentMap = TRUE}).
#' @param CS_Point.File Sample locations. Only needed when \code{jl.inputs} is
#'   not provided. See \code{\link{jl.prep}}.
#' @param pairs_to_include 1/0 vector selecting pairs. Only needed when
#'   \code{jl.inputs} is not provided.
#' @param pop2ind Per-population sample-size vector. Only needed when
#'   \code{jl.inputs} is not provided.
#' @param parallel Logical. Run Circuitscape in parallel? Default = \code{FALSE}.
#' @param cores Number of cores for parallel runs.
#' @param cholmod Logical. Use CHOLMOD solver? Default = \code{TRUE}.
#' @param precision Logical. Use single precision? Default = \code{FALSE}.
#' @param JULIA_HOME Path to the Julia \code{bin/} directory. Only needed when
#'   \code{jl.inputs} is not provided.
#' @param rm.files Logical. Remove temporary files after the run?
#'   Default = \code{TRUE}.
#' @param scratch Scratch directory for temporary files.
#' @param is_resistance Logical. Is \code{r} a resistance (vs conductance)
#'   surface? Default = \code{TRUE}.
#'
#' @return Lower-half pairwise resistance vector, a full square matrix
#'   (\code{full.mat = TRUE}), or a \code{SpatRaster} current map
#'   (\code{output = "raster"}).
#'
#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r         <- rast(nrows = 100, ncols = 100, vals = runif(1e4))
#' cs.result <- Run_CS.jl(jl.inputs = jl.inputs, r = r)
#' }
Run_CS.jl <- function(jl.inputs        = NULL,
                      r,
                      CurrentMap       = FALSE,
                      full.mat         = FALSE,
                      EXPORT.dir       = NULL,
                      output           = "matrix",
                      CS_Point.File    = NULL,
                      pairs_to_include = NULL,
                      pop2ind          = NULL,
                      parallel         = FALSE,
                      cores            = NULL,
                      cholmod          = TRUE,
                      precision        = FALSE,
                      JULIA_HOME       = NULL,
                      rm.files         = TRUE,
                      scratch          = NULL,
                      is_resistance    = TRUE) {

  wd <- getwd()
  on.exit(if (getwd() != wd) setwd(wd), add = TRUE)

  # Merge jl.inputs into local variables -------------------------------------
  if (!is.null(jl.inputs)) {
    if (!is.null(jl.inputs$rm.files)) rm.files  <- jl.inputs$rm.files
    if (!is.null(jl.inputs$scratch))  scratch   <- jl.inputs$scratch
    JULIA_HOME <- jl.inputs$JULIA_HOME

    write.files    <- jl.inputs$write.files
    write.criteria <- jl.inputs$write.criteria
  } else {
    write.files    <- NULL
    write.criteria <- NULL
  }

  if (is.null(jl.inputs) && is.null(JULIA_HOME)) {
    stop("Specify either `jl.inputs` or `JULIA_HOME` and `CS_Point.File`.")
  }
  if ((is.null(CS_Point.File) || is.null(JULIA_HOME)) && is.null(jl.inputs)) {
    stop("Both `JULIA_HOME` and `CS_Point.File` must be specified when `jl.inputs` is NULL.")
  }

  if (is.null(jl.inputs)) {
    jl.inputs <- jl.prep(
      n.Pops           = nrow(CS_Point.File),
      CS_Point.File    = CS_Point.File,
      pairs_to_include = pairs_to_include,
      JULIA_HOME       = JULIA_HOME,
      scratch          = scratch,
      parallel         = parallel,
      cores            = cores,
      cholmod          = cholmod,
      precision        = precision,
      run_test         = FALSE
    )
  }

  if (full.mat && !is.null(jl.inputs$pairs_to_include)) {
    warning(
      "The full matrix will contain -1 for excluded pairs."
    )
  }

  # Set up Julia --------------------------------------------------------------
  if (!is.null(JULIA_HOME)) {
    Sys.setenv(JULIA_BINDIR = normalizePath(JULIA_HOME))
  }
  JuliaConnectoR::juliaEval("using Circuitscape")

  # Load / validate raster ----------------------------------------------------
  if (!inherits(r, "SpatRaster")) {
    File.name <- sub("\\.(asc|tif|tiff)$", "", basename(r))
    R         <- terra::rast(r)
    names(R)  <- File.name
    asc.dir   <- r
  } else {
    R       <- r
    asc.dir <- NULL
  }

  # Export / temp directories -------------------------------------------------
  if (is.null(EXPORT.dir)) {
    if (Sys.info()[["sysname"]] == "Windows") {
      EXPORT.dir <- paste0(normalizePath(tempdir(), winslash = "/"), "/")
    } else {
      EXPORT.dir <- paste0(normalizePath(tempdir()), "/")
    }
    if (!is.null(scratch)) {
      EXPORT.dir <- paste0(normalizePath(scratch), "/")
    }
  }

  temp_rast <- rm.rast <- tempfile(pattern = "raster_",
                                   tmpdir  = tempdir(),
                                   fileext = ".asc")
  if (!is.null(scratch)) {
    temp_rast <- rm.rast <- normalizePath(
      paste0(normalizePath(scratch), "/", basename(rm.rast))
    )
  }

  tmp.name <- sub("\\.asc$", "", basename(temp_rast))

  # Current-map settings
  if (isTRUE(CurrentMap)) {
    tmp.name    <- names(R)
    MAP         <- "write_cum_cur_map_only = True"
    CURRENT.MAP <- "write_cur_maps = True"
  } else {
    MAP         <- "write_cum_cur_map_only = False"
    CURRENT.MAP <- "write_cur_maps = 0"
  }
  File.name <- names(R)

  # Prepare resistance raster -------------------------------------------------
  mx.val <- terra::global(R, "max", na.rm = TRUE)[[1]]
  if (mx.val > 1e6) R <- SCALE(R, 1, 1e6)
  R <- terra::classify(R, matrix(c(-Inf, 0, 1), ncol = 3))

  if (is.null(asc.dir)) {
    terra::writeRaster(x = R, filename = temp_rast, overwrite = TRUE)
  } else {
    temp_rast <- asc.dir
    tmp.name  <- File.name
  }

  # Write Circuitscape .ini file ----------------------------------------------
  connect <- if (jl.inputs$Neighbor.Connect == 4L) "True" else "False"

  PAIRS_TO_INCLUDE <- paste0(
    "included_pairs_file = (Browse for a file with pairs to include or exclude)"
  )
  PAIRS <- "use_included_pairs = False"

  cs.pt2 <- normalizePath(jl.inputs$CS_Point.File)
  suppressWarnings({
    OUT   <- paste0("output_file = ",
                    normalizePath(paste0(EXPORT.dir, tmp.name, ".out")))
    BATCH <- normalizePath(paste0(EXPORT.dir, tmp.name, ".ini"))
  })

  write.CS_4.0(
    BATCH            = BATCH,
    OUT              = OUT,
    HABITAT          = paste0("habitat_file = ", temp_rast),
    LOCATION.FILE    = paste0("point_file = ", cs.pt2),
    CONNECTION       = paste0("connect_four_neighbors_only =", connect),
    MAP              = MAP,
    CURRENT.MAP      = CURRENT.MAP,
    PAIRS_TO_INCLUDE = PAIRS_TO_INCLUDE,
    PAIRS            = PAIRS,
    PARALLELIZE      = jl.inputs$parallel,
    CORES            = jl.inputs$cores,
    solver           = jl.inputs$solver,
    precision        = jl.inputs$precision,
    silent           = jl.inputs$silent,
    is_resistance    = is_resistance
  )

  # Execute Circuitscape.jl ---------------------------------------------------
  rt <- NULL

  if (!is.null(write.criteria)) {
    t1  <- proc.time()[3]
    out <- JuliaConnectoR::juliaCall(
      "compute",
      normalizePath(paste0(EXPORT.dir, tmp.name, ".ini"))
    )[-1, -1]
    rt  <- proc.time()[3] - t1
  } else {
    out <- JuliaConnectoR::juliaCall(
      "compute",
      normalizePath(paste0(EXPORT.dir, tmp.name, ".ini"))
    )[-1, -1]
  }

  # Return raster output ------------------------------------------------------
  if (output == "raster" && isTRUE(CurrentMap)) {
    cur_file <- normalizePath(
      paste0(EXPORT.dir, "/", tmp.name, "_cum_curmap.asc")
    )
    if (!file.exists(cur_file) && !is.null(scratch)) {
      cur_file <- normalizePath(paste0(scratch, "/", tmp.name, "_cum_curmap.asc"))
    }
    rast_out      <- terra::rast(cur_file)
    names(rast_out) <- File.name

    if (isTRUE(rm.files)) {
      del <- c(
        list.files(normalizePath(EXPORT.dir), pattern = tmp.name,
                   full.names = TRUE, all.files = TRUE),
        list.files(normalizePath(EXPORT.dir), pattern = basename(temp_rast),
                   full.names = TRUE, all.files = TRUE)
      )
      invisible(sapply(del, unlink, force = TRUE))
    }
    return(rast_out)
  }

  # Return distance matrix or vector ------------------------------------------
  if (!full.mat) {
    if (!is.null(jl.inputs$pairs_to_include)) {
      cs.matrix <- lower(out)
      cs.matrix <- cs.matrix[jl.inputs$keep != 0]
    } else if (!is.null(pop2ind)) {
      cs.matrix <- expand.mat(out, pop2ind)
    } else {
      cs.matrix <- lower(out)
      cs.matrix <- cs.matrix[cs.matrix != -1]
    }
  } else {
    if (!is.null(jl.inputs$pairs_to_include)) {
      mat <- matrix(0, nrow = jl.inputs$n.Pops, ncol = jl.inputs$n.Pops)
      mat[lower.tri(mat)] <- jl.inputs$keep
      mat  <- mat + t(mat) - diag(diag(mat))
      cs.matrix <- mat * out
      cs.matrix <- ifelse(cs.matrix == 0, -1, cs.matrix)
      diag(cs.matrix) <- 0
    } else if (!is.null(pop2ind)) {
      cs.matrix <- expand.mat(out, pop2ind, format = "matrix")
    } else {
      cs.matrix <- out
    }
  }

  if (any(cs.matrix == 0, na.rm = TRUE)) {
    message(
      "Zero values generated by Circuitscape.\n",
      "Check whether multiple sample points fall in the same raster cell."
    )
  }
  if (any(is.na(cs.matrix))) {
    cs.matrix[is.na(cs.matrix)] <- 0
    message(
      "NA values generated by Circuitscape.\n",
      "Check whether multiple sample points fall in the same raster cell."
    )
  }

  # Optionally copy .ini / .asc files ----------------------------------------
  if (!is.null(write.files)) {
    do_write <- is.null(write.criteria) ||
      (!is.null(rt) && rt > write.criteria)
    if (do_write) {
      file.copy(
        c(normalizePath(paste0(EXPORT.dir, tmp.name, ".ini")),
          normalizePath(temp_rast)),
        normalizePath(write.files)
      )
    }
  }

  # Remove temporary files ----------------------------------------------------
  if (isTRUE(rm.files)) {
    del <- list.files(normalizePath(EXPORT.dir),
                      pattern   = tmp.name,
                      full.names = TRUE,
                      all.files  = TRUE)
    invisible(sapply(del, unlink, force = TRUE))

    if (!is.null(scratch)) {
      del2 <- list.files(normalizePath(scratch),
                         pattern   = tmp.name,
                         full.names = TRUE,
                         all.files  = TRUE)
      invisible(sapply(del2, unlink, force = TRUE))
      unlink(rm.rast, force = TRUE)
    }
  }

  # Clean up stale .asc files in EXPORT.dir -----------------------------------
  asc.files  <- list.files(EXPORT.dir, pattern = "\\.asc$",
                            full.names = TRUE, all.files = TRUE)
  if (length(asc.files) > 0) {
    time.diff <- difftime(Sys.time(), file.mtime(asc.files), units = "secs")
    unlink(asc.files[time.diff > 30], force = TRUE)
  }

  return(cs.matrix)
}
