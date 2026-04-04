#' Prepare inputs for running Circuitscape via Julia
#'
#' Creates all inputs required by the optimization functions when using
#' Circuitscape run through Julia (via \pkg{JuliaConnectoR}).
#'
#' @param n.Pops The number of populations being assessed.
#' @param response Vector of pairwise genetic distances (lower half of pairwise
#'   matrix). If using \code{pairs_to_include}, provide the full lower-half
#'   vector; it will be trimmed automatically.
#' @param CS_Point.File A \code{\link[terra]{SpatVector}} of sample point
#'   locations.
#' @param covariates Data frame of additional covariates for the MLPE model.
#' @param formula R formula for the fixed effects of the MLPE model (e.g.,
#'   \code{response ~ covariate}). The \code{response} term uses the values in
#'   \code{response}; covariate names must match \code{covariates} columns.
#' @param JULIA_HOME Path to the folder containing the Julia binary. If
#'   \code{NULL} Julia must already be on the system PATH.
#' @param Neighbor.Connect Connectivity scheme for Circuitscape: 4 or 8
#'   (Default = 8).
#' @param pairs_to_include Integer vector of 1 (keep) or 0 (drop) for each
#'   pairwise observation. Default = \code{NULL} (use all pairs).
#' @param pop2ind Integer vector giving the number of sampled individuals per
#'   population. Default = \code{NULL}.
#' @param nb Maximum distance (in coordinate units) for treating locations as
#'   the same neighborhood. Default = \code{NULL}.
#' @param parallel Logical. Run Circuitscape in parallel? Default = \code{FALSE}.
#' @param cores Number of cores for parallel processing (only used when
#'   \code{parallel = TRUE}).
#' @param cholmod Logical. Use the CHOLMOD direct solver? Default = \code{TRUE}.
#'   Cannot be used with single precision (\code{precision = TRUE}).
#' @param precision Logical. Use experimental single-precision mode?
#'   Default = \code{FALSE}.
#' @param run_test Logical. Test the Julia/Circuitscape connection before
#'   continuing? Default = \code{TRUE}.
#' @param write.files Directory path. If provided, the \code{.ini} and
#'   \code{.asc} files used in each Circuitscape run are saved there. If the
#'   directory does not exist, it will be created automatically.
#' @param write.criteria Minimum run time (seconds) threshold for writing files
#'   when \code{write.files} is set. If \code{NULL}, all runs are written.
#' @param silent Logical. Suppress Circuitscape progress output?
#'   Default = \code{TRUE}.
#' @param scratch Scratch directory for temporary files when the default
#'   \code{tempdir()} is not writable.
#' @param rm.files Logical. Remove temporary files after each Julia run?
#'   Default = \code{TRUE}.
#'
#' @return A named list of inputs required by the optimization functions.
#'
#' @details
#' Julia must be installed on your system. On the first call, the
#' \code{Circuitscape.jl} package will be downloaded if it is not already
#' present (see \url{https://github.com/Circuitscape/Circuitscape.jl}).
#'
#' Set \code{JULIA_HOME} to the \code{bin/} subdirectory of your Julia
#' installation (e.g., \code{"/usr/local/julia/bin"}).
#'
#' The CHOLMOD solver (\code{cholmod = TRUE}) uses a direct Cholesky
#' decomposition that can be much faster than the default algebraic multigrid
#' solver for moderate problem sizes, but may consume large amounts of memory
#' for very large problems. It cannot be used with single precision.
#'
#' Point locations are supplied as a \code{SpatVector} and written to a
#' temporary Circuitscape point file internally.
#'
#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @examples
#' \dontrun{
#' pts <- terra::vect(sample_pops[[1]], type = "points")
#'
#' jl.inputs <- jl.prep(
#'   n.Pops = nrow(sample_pops[[1]]),
#'   response = lower(Dc_list[[1]]),
#'   CS_Point.File = pts,
#'   JULIA_HOME = Sys.getenv("JULIA_BINDIR"),
#'   run_test = FALSE
#' )
#' }
jl.prep <- function(n.Pops,
                    response         = NULL,
                    CS_Point.File,
                    covariates       = NULL,
                    formula          = NULL,
                    JULIA_HOME       = NULL,
                    Neighbor.Connect = 8,
                    pairs_to_include = NULL,
                    pop2ind          = NULL,
                    nb               = NULL,
                    parallel         = FALSE,
                    cores            = NULL,
                    cholmod          = TRUE,
                    precision        = FALSE,
                    run_test         = TRUE,
                    write.files      = NULL,
                    write.criteria   = NULL,
                    silent           = TRUE,
                    scratch          = NULL,
                    rm.files         = TRUE) {

  # Input validation ----------------------------------------------------------

  if (!is.null(covariates) && !is.data.frame(covariates)) {
    stop("'covariates' must be a data frame.")
  }
  if (!is.null(covariates) && !is.null(response) &&
      nrow(covariates) != length(response)) {
    stop("'response' and 'covariates' must have the same number of observations.")
  }
  cs_points <- .validate_spatvector_points(CS_Point.File, arg = "CS_Point.File")

  if (!is.null(write.files)) {
    write.files <- paste0(
      sub("[/\\\\]+$", "", normalizePath(write.files, winslash = "/", mustWork = FALSE)),
      "/"
    )
    if (!dir.exists(write.files)) {
      dir.create(write.files, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(write.files)) {
        stop("Failed to create 'write.files' directory: ", write.files)
      }
      message("Created 'write.files' directory: ", write.files)
    }
  }

  # Julia setup ---------------------------------------------------------------
  if (!is.null(JULIA_HOME)) {
    if (!dir.exists(JULIA_HOME))
      stop("Specified JULIA_HOME directory does not exist: ", JULIA_HOME)
    Sys.setenv(JULIA_BINDIR = normalizePath(JULIA_HOME))
  }
  JULIA_HOME <- normalizePath(JULIA_HOME %||% Sys.getenv("JULIA_BINDIR"))

  # Check Circuitscape.jl is installed
  cs_installed <- tryCatch({
    JuliaConnectoR::juliaEval("using Circuitscape")
    TRUE
  }, error = function(e) FALSE)

  if (!cs_installed) {
    stop(
      "The Julia Circuitscape package is not installed.\n",
      "Install it from: https://github.com/Circuitscape/Circuitscape.jl"
    )
  }

  # Solver and precision ------------------------------------------------------
  if (isTRUE(precision)) {
    precision <- "single"
  } else {
    precision <- "None"
  }

  if (isTRUE(cholmod) && precision == "single") {
    stop(
      "CHOLMOD solver requires double precision. ",
      "Set either `cholmod = FALSE` or `precision = FALSE`."
    )
  }

  solver <- if (isTRUE(cholmod)) "cholmod" else NULL

  # Determine temp directory --------------------------------------------------
  if (Sys.info()[["sysname"]] == "Windows") {
    td <- paste0(normalizePath(tempdir(), winslash = "/"), "/")
  } else {
    td <- paste0(normalizePath(tempdir()), "/")
  }
  if (!is.null(scratch)) {
    td <- paste0(normalizePath(scratch), "/")
  }

  # Parse CS_Point.File -------------------------------------------------------
  cs_coords <- .point_coords(cs_points, arg = "CS_Point.File")
  if (n.Pops != nrow(cs_coords)) {
    stop("'n.Pops' (", n.Pops, ") does not equal the number of sample locations (",
         nrow(cs_coords), ").")
  }
  site      <- seq_len(nrow(cs_points))
  cs.txt    <- data.frame(site, cs_coords)
  write.table(cs.txt,
              file      = paste0(td, "sample_pts.txt"),
              col.names = FALSE,
              row.names = FALSE)
  CS_Point.File <- paste0(td, "sample_pts.txt")

  # Test run ------------------------------------------------------------------
  if (isTRUE(run_test)) {
    message("Testing Circuitscape via Julia...")

    test_pts <- ResistanceGA2::samples[1:5, ]
    write.table(test_pts,
                file      = paste0(td, "samples.txt"),
                quote     = FALSE,
                sep       = "\t",
                row.names = FALSE,
                col.names = FALSE)

    temp.ini  <- tempfile(pattern = "", tmpdir = td, fileext = ".ini")
    tmp.name  <- sub("\\.ini$", "", basename(temp.ini))

    terra::writeRaster(
      x         = ResistanceGA2::resistance_surfaces[["continuous"]],
      filename  = paste0(td, tmp.name, ".asc"),
      overwrite = TRUE
    )

    write.CS_4.0(
      BATCH          = temp.ini,
      OUT            = paste0("output_file = ", td, tmp.name, ".out"),
      HABITAT        = paste0("habitat_file = ", td, tmp.name, ".asc"),
      LOCATION.FILE  = paste0("point_file = ", td, "samples.txt"),
      PARALLELIZE    = FALSE,
      CORES          = NULL,
      solver         = "cholmod",
      precision      = "None",
      silent         = silent
    )

    wd  <- getwd()
    invisible(capture.output(
      out <- JuliaConnectoR::juliaCall("compute", normalizePath(temp.ini))[-1, -1],
      type = "message"
    ))
    if (getwd() != wd) setwd(wd)

    if (nrow(out) == 5L) {
      message("Test passed.")
    } else {
      stop("Circuitscape Julia test failed.")
    }

    if (isTRUE(rm.files)) {
      del <- list.files(td, pattern = tmp.name, full.names = TRUE)
      invisible(sapply(del, unlink, force = TRUE))
    }
  }

  # ID / ZZ matrices ----------------------------------------------------------
  if (!is.null(nb)) {
    ID <- To.From.ID(sampled_pops = n.Pops,
                     pop_n        = pop2ind,
                     spLoc        = cs_points,
                     nb           = nb)
  } else {
    ID <- To.From.ID(sampled_pops = n.Pops,
                     pop_n        = pop2ind,
                     spLoc        = NULL,
                     nb           = nb)
  }

  suppressWarnings(ZZ <- ZZ.mat(ID))

  # Build data frame and formula ----------------------------------------------
  fmla <- formula
  df   <- NULL

  if (!is.null(response)) {
    if (!is.null(covariates)) {
      if (!is.null(pop2ind)) {
        keep.vec  <- expand.keep(pop2ind)
        gd.mat    <- matrix(0, n.Pops, n.Pops)
        cov.list  <- lapply(seq_len(ncol(covariates)), function(i) {
          mat <- gd.mat
          mat[lower.tri(mat)] <- covariates[, i]
          expand.mat(mat, pop2ind)
        })
        names(cov.list) <- names(covariates)
        cov.df    <- as.data.frame(cov.list)
        covariates <- cov.df
        response  <- response[keep.vec == 1]

        if (!is.null(nb)) {
          df <- data.frame(gd = response, cov.df,
                           pop = ID$pop, grp = ID$pop1.pop,
                           cor.grp = ID$cor.grp)
        } else {
          df <- data.frame(gd = response, cov.df,
                           pop = ID$pop, grp = ID$pop1.pop)
        }
      } else {
        if (!is.null(nb)) {
          df <- data.frame(gd = response, covariates,
                           pop = ID$pop1, cor.grp = ID$cor.grp)
        } else {
          df <- data.frame(gd = response, covariates, pop = ID$pop1)
        }
      }
    } else {
      if (!is.null(pop2ind)) {
        keep.vec <- expand.keep(pop2ind)
        response <- response[keep.vec == 1]
        if (!is.null(nb)) {
          df <- data.frame(gd = response, pop = ID$pop,
                           grp = ID$pop1.pop, cor.grp = ID$cor.grp)
        } else {
          df <- data.frame(gd = response, pop = ID$pop, grp = ID$pop1.pop)
        }
      } else {
        if (!is.null(nb)) {
          df <- data.frame(gd = response, pop = ID$pop1, cor.grp = ID$cor.grp)
        } else {
          df <- data.frame(gd = response, pop = ID$pop1)
        }
      }
    }
  }

  # Formula
  if (!is.null(fmla)) {
    if (!is.null(pop2ind)) {
      if (!is.null(nb)) {
        fmla <- update(fmla, gd ~ . + cd + (1 | pop) + (1 | grp) + (1 | cor.grp))
      } else {
        fmla <- update(fmla, gd ~ . + cd + (1 | pop) + (1 | grp))
      }
    } else {
      if (!is.null(nb)) {
        fmla <- update(fmla, gd ~ . + cd + (1 | pop) + (1 | cor.grp))
      } else {
        fmla <- update(fmla, gd ~ . + cd + (1 | pop))
      }
    }
  } else {
    if (!is.null(pop2ind)) {
      fmla <- if (!is.null(nb)) {
        gd ~ cd + (1 | pop) + (1 | grp) + (1 | cor.grp)
      } else {
        gd ~ cd + (1 | pop) + (1 | grp)
      }
    } else {
      fmla <- if (!is.null(nb)) {
        gd ~ cd + (1 | pop) + (1 | cor.grp)
      } else {
        gd ~ cd + (1 | pop)
      }
    }
  }

  # Pairs to include ----------------------------------------------------------
  pairs_to_include.file <- NULL
  keep                  <- pairs_to_include
  if (is.null(keep)) keep <- rep(1L, nrow(ID))

  ID.keep          <- ID
  covariates.keep  <- covariates
  ZZ.keep          <- ZZ
  response.keep    <- response

  if (!is.null(pairs_to_include)) {
    keep    <- pairs_to_include
    df      <- df[keep == 1, , drop = FALSE]
    ZZ.keep <- ZZ[, keep == 1, drop = FALSE]
    pop     <- ID$pop1[keep == 1]
    ZZ.keep <- ZZ.keep[rownames(ZZ.keep) %in% unique(as.character(pop)), , drop = FALSE]

    response.keep   <- response[keep == 1]
    covariates.keep <- if (!is.null(covariates)) covariates[keep == 1, , drop = FALSE] else NULL
    ID.keep         <- ID[keep == 1, , drop = FALSE]

    keep.id <- ID[keep == 1, , drop = FALSE]
    names(keep.id) <- c("mode", "include")
    write.table(keep.id,
                file      = paste0(td, "include_pairs.txt"),
                col.names = TRUE,
                row.names = FALSE,
                quote     = FALSE)
    pairs_to_include.file <- paste0(td, "include_pairs.txt")
  }

  # Check slope ---------------------------------------------------------------
  if (!is.null(response) && !is.null(cs_coords)) {
    cd  <- c(dist(cs_coords))[keep == 1]
    m   <- lm(gd ~ cd, data = df)
    if (coef(m)[2] < 0) {
      warning(
        "Genetic distance decreases with Euclidean distance. ",
        "This is likely to result in a failed optimization.\n",
        "Consider subtracting your values from 1 to reverse the relationship."
      )
    }
  }

  if (!is.null(df)) {
    df <- .mlpe_attach_workflow_pairs(df, ID.keep)
  }

  # Return --------------------------------------------------------------------
  list(
    ID                    = ID.keep,
    ZZ                    = ZZ.keep,
    df                    = df,
    response              = response.keep,
    covariates            = covariates.keep,
    formula               = fmla,
    CS_Point.File         = CS_Point.File,
    Neighbor.Connect      = Neighbor.Connect,
    n.Pops                = n.Pops,
    pairs_to_include      = pairs_to_include.file,
    pop2ind               = pop2ind,
    nb                    = nb,
    parallel              = parallel,
    cores                 = cores,
    solver                = solver,
    precision             = precision,
    JULIA_HOME            = JULIA_HOME,
    write.files           = write.files,
    write.criteria        = write.criteria,
    silent                = silent,
    scratch               = scratch,
    keep                  = keep,
    response.all          = response,
    ID.all                = ID,
    ZZ.all                = ZZ,
    covariates.all        = covariates,
    rm.files              = rm.files,
    opt.input             = "jl"
  )
}

# Null-coalescing operator (unexported helper)
`%||%` <- function(a, b) if (!is.null(a)) a else b
