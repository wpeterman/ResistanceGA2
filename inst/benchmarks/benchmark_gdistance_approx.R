# Compare exact gdistance commute runs with aggregate approximations.
#
# Example:
# Rscript inst/benchmarks/benchmark_gdistance_approx.R
#
# Optional environment variables:
# RGA_APPROX_SIZES="1000"
# RGA_APPROX_NPTS="100"
# RGA_APPROX_REPS="10"
# RGA_APPROX_FACTORS="2,4"

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

make_case <- function(size, npts, seed, directions = 8) {
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

  list(
    raster = r,
    exact = ResistanceGA2::gdist.prep(
      n.Pops = npts,
      samples = pts,
      method = "commuteDistance",
      directions = directions
    )
  )
}

time_call <- function(expr) {
  gc()
  elapsed <- system.time(value <- force(expr))[["elapsed"]]
  list(value = value, elapsed = elapsed)
}

benchmark_one <- function(size, npts, rep, factors, directions = 8) {
  case <- make_case(
    size = size,
    npts = npts,
    seed = (size * 1000) + (npts * 10) + rep,
    directions = directions
  )

  exact <- time_call(
    ResistanceGA2::Run_gdistance(case$exact, case$raster)
  )
  exact_vec <- as.numeric(exact$value)
  rows <- list(data.frame(
    engine = "exact",
    raster_cells = size * size,
    raster_nrow = size,
    n_points = npts,
    directions = directions,
    aggregate_factor = 1L,
    replicate = rep,
    elapsed_sec = exact$elapsed,
    output_cor = 1,
    output_spearman = 1,
    median_ratio = 1,
    rel_mae = 0,
    stringsAsFactors = FALSE
  ))

  for (factor in factors) {
    approx_inputs <- case$exact
    approx_inputs$commute.approx <- "aggregate"
    approx_inputs$approx.factor <- as.integer(factor)
    approx_inputs$approx.scale <- TRUE

    approx <- time_call(
      ResistanceGA2::Run_gdistance(approx_inputs, case$raster)
    )
    approx_vec <- as.numeric(approx$value)
    rows[[length(rows) + 1L]] <- data.frame(
      engine = paste0("aggregate", factor),
      raster_cells = size * size,
      raster_nrow = size,
      n_points = npts,
      directions = directions,
      aggregate_factor = factor,
      replicate = rep,
      elapsed_sec = approx$elapsed,
      output_cor = suppressWarnings(stats::cor(exact_vec, approx_vec)),
      output_spearman = suppressWarnings(stats::cor(
        exact_vec,
        approx_vec,
        method = "spearman"
      )),
      median_ratio = stats::median(approx_vec) / stats::median(exact_vec),
      rel_mae = mean(abs(approx_vec - exact_vec)) / mean(exact_vec),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

main <- function() {
  sizes <- parse_int_env("RGA_APPROX_SIZES", 1000L)
  npts <- parse_int_env("RGA_APPROX_NPTS", 100L)[[1]]
  reps <- parse_int_env("RGA_APPROX_REPS", 10L)[[1]]
  factors <- parse_int_env("RGA_APPROX_FACTORS", c(2L, 4L))
  directions <- parse_int_env("RGA_APPROX_DIRECTIONS", 8L)[[1]]

  message("Benchmark sizes: ", paste(sizes, collapse = ", "))
  message("Points: ", npts, "; reps: ", reps, "; directions: ", directions)
  message("Aggregate factors: ", paste(factors, collapse = ", "))

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
            factors = factors,
            directions = directions
          )
        })
      )
    })
  )
  rownames(results) <- NULL

  print(results)

  aggregate_summary <- aggregate(
    cbind(elapsed_sec, output_cor, output_spearman, median_ratio, rel_mae) ~
      raster_nrow + n_points + engine,
    data = results,
    FUN = function(x) c(mean = mean(x), median = stats::median(x),
                        min = min(x), max = max(x))
  )
  aggregate_summary <- do.call(data.frame, aggregate_summary)
  names(aggregate_summary) <- sub("\\.", "_", names(aggregate_summary))

  cat("\nSummary:\n")
  print(aggregate_summary)

  invisible(results)
}

if (identical(environment(), globalenv())) {
  main()
}
