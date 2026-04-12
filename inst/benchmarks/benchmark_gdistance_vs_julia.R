# Compare optimized gdistance commute runs with Julia/Circuitscape.
#
# Example:
# Rscript inst/benchmarks/benchmark_gdistance_vs_julia.R
#
# Optional environment variables:
# RGA_BENCH_SIZES="50,100,250"
# RGA_BENCH_NPTS="20"
# RGA_BENCH_REPS="1"
# RGA_BENCH_WARMUP="TRUE"
# RGA_BENCH_KEEP_PAIRS="300"
# JULIA_BINDIR="C:/path/to/Julia/bin"

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") &&
      grepl("^Package:\\s+ResistanceGA2\\s*$", readLines("DESCRIPTION", n = 1))) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      stop("Install devtools or run this benchmark from an installed ResistanceGA2 package.")
    }
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(ResistanceGA2)
  }
  library(terra)
})

parse_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

parse_logical_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  toupper(value) %in% c("TRUE", "T", "1", "YES", "Y")
}

find_julia_bindir <- function() {
  bindir <- Sys.getenv("JULIA_BINDIR", unset = "")
  if (nzchar(bindir) && file.exists(file.path(bindir, "julia.exe"))) {
    return(normalizePath(bindir, winslash = "/", mustWork = TRUE))
  }

  candidates <- Sys.glob(file.path(
    Sys.getenv("LOCALAPPDATA"),
    "Programs",
    "Julia-*",
    "bin",
    "julia.exe"
  ))
  if (length(candidates)) {
    versions <- sub("^Julia-", "", basename(dirname(dirname(candidates))))
    newest <- candidates[[order(numeric_version(versions), decreasing = TRUE)[[1]]]]
    return(normalizePath(dirname(newest), winslash = "/", mustWork = TRUE))
  }

  julia <- Sys.which("julia")
  if (nzchar(julia) && !grepl("WindowsApps", julia, fixed = TRUE)) {
    return(normalizePath(dirname(julia), winslash = "/", mustWork = TRUE))
  }

  stop("Julia not found. Set JULIA_BINDIR to a Julia bin directory.")
}

make_case <- function(size, npts, seed, directions = 8, keep_pairs = NA_integer_) {
  r <- terra::rast(
    nrows = size, ncols = size,
    xmin = 0, xmax = size,
    ymin = 0, ymax = size
  )
  set.seed(seed)
  terra::values(r) <- runif(terra::ncell(r), 1, 100)

  inner_cells <- expand.grid(
    col = seq.int(2L, size - 1L),
    row = seq.int(2L, size - 1L)
  )
  if (npts > nrow(inner_cells)) {
    stop("'npts' is larger than the number of available inner raster cells.")
  }
  sampled_cells <- inner_cells[sample.int(nrow(inner_cells), npts), , drop = FALSE]
  coords <- cbind(
    sampled_cells$col - 0.5,
    size - sampled_cells$row + 0.5
  )
  pts <- terra::vect(coords, type = "points")
  n_pairs <- npts * (npts - 1L) / 2L
  keep <- NULL
  if (!is.na(keep_pairs)) {
    keep_pairs <- min(as.integer(keep_pairs), n_pairs)
    keep <- integer(n_pairs)
    keep[sample.int(n_pairs, keep_pairs)] <- 1L
  }

  list(
    raster = r,
    points = pts,
    gdist = ResistanceGA2::gdist.prep(
      n.Pops = npts,
      samples = pts,
      method = "commuteDistance",
      directions = directions,
      keep = keep
    ),
    keep = keep
  )
}

time_call <- function(expr) {
  gc()
  elapsed <- system.time(value <- force(expr))[["elapsed"]]
  list(value = value, elapsed = elapsed)
}

vectorize_gdistance <- function(x) {
  as.numeric(x)
}

benchmark_one <- function(size, npts, rep, julia_home, scratch, directions = 8,
                          keep_pairs = NA_integer_) {
  case <- make_case(
    size = size,
    npts = npts,
    seed = (size * 1000) + (npts * 10) + rep,
    directions = directions,
    keep_pairs = keep_pairs
  )

  jl <- ResistanceGA2::jl.prep(
    n.Pops = npts,
    CS_Point.File = case$points,
    JULIA_HOME = julia_home,
    Neighbor.Connect = directions,
    pairs_to_include = case$keep,
    run_test = FALSE,
    scratch = scratch,
    silent = TRUE
  )

  gdist <- time_call(
    ResistanceGA2::Run_gdistance(case$gdist, case$raster)
  )
  julia <- time_call(
    ResistanceGA2::Run_CS.jl(jl.inputs = jl, r = case$raster, full.mat = FALSE)
  )

  gdist_vec <- vectorize_gdistance(gdist$value)
  julia_vec <- as.numeric(julia$value)

  data.frame(
    engine = c("gdistance_fast", "julia_circuitscape"),
    raster_cells = size * size,
    raster_nrow = size,
    raster_ncol = size,
    n_points = npts,
    directions = directions,
    kept_pairs = if (is.null(case$keep)) length(gdist_vec) else sum(case$keep),
    replicate = rep,
    elapsed_sec = c(gdist$elapsed, julia$elapsed),
    output_length = c(length(gdist_vec), length(julia_vec)),
    output_cor = suppressWarnings(stats::cor(gdist_vec, julia_vec)),
    output_spearman = suppressWarnings(stats::cor(gdist_vec, julia_vec, method = "spearman")),
    median_value = c(stats::median(gdist_vec), stats::median(julia_vec)),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  sizes <- parse_int_env("RGA_BENCH_SIZES", c(50L, 100L, 250L))
  npts <- parse_int_env("RGA_BENCH_NPTS", 20L)[[1]]
  reps <- parse_int_env("RGA_BENCH_REPS", 1L)[[1]]
  directions <- parse_int_env("RGA_BENCH_DIRECTIONS", 8L)[[1]]
  warmup <- parse_logical_env("RGA_BENCH_WARMUP", TRUE)
  keep_pairs <- parse_int_env("RGA_BENCH_KEEP_PAIRS", NA_integer_)[[1]]
  julia_home <- find_julia_bindir()
  scratch <- file.path(tempdir(), "rga-julia-benchmark")
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)

  message("Julia bin dir: ", julia_home)
  message("Benchmark sizes: ", paste(sizes, collapse = ", "))
  message("Points: ", npts, "; reps: ", reps, "; directions: ", directions)
  message("Kept pairs: ", ifelse(is.na(keep_pairs), "all", keep_pairs))
  message("Warmup: ", warmup)

  if (isTRUE(warmup)) {
    invisible(benchmark_one(
      size = 15L,
      npts = min(npts, 5L),
      rep = 0L,
      julia_home = julia_home,
      scratch = scratch,
      directions = directions,
      keep_pairs = if (is.na(keep_pairs)) NA_integer_ else min(keep_pairs, 10L)
    ))
  }

  results <- do.call(
    rbind,
    lapply(sizes, function(size) {
      do.call(
        rbind,
        lapply(seq_len(reps), function(rep) {
          benchmark_one(
            size = size,
            npts = npts,
            rep = rep,
            julia_home = julia_home,
            scratch = scratch,
            directions = directions,
            keep_pairs = keep_pairs
          )
        })
      )
    })
  )

  results <- results[order(results$raster_cells, results$replicate, results$engine), ]
  rownames(results) <- NULL

  print(results)

  wide <- reshape(
    results[, c("raster_nrow", "n_points", "kept_pairs", "replicate", "engine", "elapsed_sec")],
    idvar = c("raster_nrow", "n_points", "kept_pairs", "replicate"),
    timevar = "engine",
    direction = "wide"
  )
  wide$speedup_vs_julia <- wide$elapsed_sec.julia_circuitscape /
    wide$elapsed_sec.gdistance_fast

  cat("\nSpeedup summary (Julia / gdistance elapsed):\n")
  print(wide)

  aggregate_summary <- aggregate(
    elapsed_sec ~ raster_nrow + n_points + kept_pairs + engine,
    data = results,
    FUN = function(x) c(mean = mean(x), median = stats::median(x), min = min(x), max = max(x))
  )
  aggregate_summary <- do.call(
    data.frame,
    aggregate_summary
  )
  names(aggregate_summary) <- sub("^elapsed_sec\\.", "", names(aggregate_summary))

  cat("\nElapsed-time summary by engine:\n")
  print(aggregate_summary)

  invisible(results)
}

if (identical(environment(), globalenv())) {
  main()
}
