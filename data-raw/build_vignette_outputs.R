args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 0) {
  stop("Run this script with Rscript so the project root can be determined.")
}

script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
project_root <- dirname(dirname(script_path))
setwd(project_root)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The 'pkgload' package is required to build vignette outputs.")
}

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
    stop("No Julia executable found. Set JULIA_BINDIR or add Julia to PATH.")
  }

  normalizePath(dirname(candidates[[1]]), winslash = "/", mustWork = TRUE)
}

pkgload::load_all(".", quiet = TRUE)

pts <- terra::vect(samples[, 2:3], type = "points")
surface <- resistance_surfaces[["continuous"]]

message("Building gdistance vignette outputs...")
gdist_inputs <- gdist.prep(
  n.Pops = nrow(samples),
  samples = pts,
  method = "costDistance"
)

gdistance_commute <- Run_gdistance(
  gdist.inputs = gdist_inputs,
  r = surface
)
if (inherits(gdistance_commute, "dist")) {
  gdistance_commute <- as.matrix(gdistance_commute)
}

message("Building Circuitscape vignette outputs...")
bindir <- find_julia_bindir()
Sys.setenv(JULIA_BINDIR = bindir)

tryCatch(
  JuliaConnectoR::juliaEval("using Circuitscape"),
  error = function(e) {
    stop("Circuitscape.jl is not available in Julia: ", conditionMessage(e))
  }
)

jl_inputs <- jl.prep(
  n.Pops = nrow(samples),
  CS_Point.File = pts,
  JULIA_HOME = bindir,
  run_test = FALSE,
  silent = TRUE,
  rm.files = FALSE
)

circuitscape_matrix <- Run_CS.jl(
  jl.inputs = jl_inputs,
  r = surface,
  CurrentMap = FALSE,
  output = "matrix"
)

export_dir <- tempfile("rga2-vignette-current-")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(export_dir, recursive = TRUE, force = TRUE), add = TRUE)

vignette_current_map <- Run_CS.jl(
  jl.inputs = jl_inputs,
  r = surface,
  CurrentMap = TRUE,
  output = "raster",
  EXPORT.dir = export_dir,
  rm.files = FALSE
)

vignette_pairwise_outputs <- list(
  gdistance_commute = gdistance_commute,
  circuitscape_matrix = circuitscape_matrix
)

save(
  vignette_pairwise_outputs,
  file = file.path("data", "vignette_pairwise_outputs.rda"),
  compress = "bzip2"
)

vignette_current_map <- terra::wrap(vignette_current_map)

save(
  vignette_current_map,
  file = file.path("data", "vignette_current_map.rda"),
  compress = "bzip2"
)

message("Saved data/vignette_pairwise_outputs.rda and data/vignette_current_map.rda")
