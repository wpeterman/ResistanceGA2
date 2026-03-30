find_julia_bindir <- function() {
  bindir <- Sys.getenv("JULIA_BINDIR", unset = "")
  if (nzchar(bindir) && dir.exists(bindir)) {
    return(normalizePath(bindir, winslash = "/", mustWork = TRUE))
  }

  candidates <- character()
  julia <- Sys.which("julia")
  if (nzchar(julia)) {
    candidates <- c(candidates, julia)
  }

  if (.Platform$OS.type == "windows") {
    where_hits <- tryCatch(
      system2("where.exe", "julia", stdout = TRUE, stderr = FALSE),
      error = function(e) character()
    )
    candidates <- c(candidates, where_hits)
  }

  candidates <- unique(candidates[nzchar(candidates)])
  candidates <- candidates[file.exists(candidates)]

  if (.Platform$OS.type == "windows") {
    non_alias <- candidates[!grepl("WindowsApps", candidates, ignore.case = TRUE)]
    if (length(non_alias) > 0) {
      candidates <- c(non_alias, setdiff(candidates, non_alias))
    }
  }

  if (length(candidates) == 0) {
    return("")
  }

  normalizePath(dirname(candidates[[1]]), winslash = "/", mustWork = TRUE)
}

test_that("Run_CS.jl auto-creates EXPORT.dir and returns a current map", {
  testthat::skip_if_not_installed("JuliaConnectoR")

  bindir <- find_julia_bindir()
  testthat::skip_if(bindir == "", "Julia executable not available.")

  Sys.setenv(JULIA_BINDIR = bindir)
  cs_ready <- tryCatch({
    JuliaConnectoR::juliaEval("using Circuitscape")
    TRUE
  }, error = function(e) FALSE)
  testthat::skip_if_not(cs_ready, "Circuitscape.jl is not available in Julia.")

  pts <- terra::vect(sample_pops[[1]][1:5, ], type = "points")
  gd <- lower(Dc_list[[1]][1:5, 1:5])

  jl_inputs <- jl.prep(
    n.Pops = 5,
    response = gd,
    CS_Point.File = pts,
    JULIA_HOME = bindir,
    run_test = FALSE,
    silent = TRUE
  )

  export_dir <- tempfile("rga2-jl-current-")
  on.exit(unlink(export_dir, recursive = TRUE, force = TRUE), add = TRUE)

  out <- suppressMessages(
    Run_CS.jl(
      jl.inputs = jl_inputs,
      r = raster_orig[["cont_orig"]],
      CurrentMap = TRUE,
      output = "raster",
      EXPORT.dir = export_dir
    )
  )

  expect_true(dir.exists(export_dir))
  expect_true(inherits(out, "SpatRaster"))
  expect_equal(terra::nlyr(out), 1L)
  expect_identical(names(out), "cont_orig")
})
