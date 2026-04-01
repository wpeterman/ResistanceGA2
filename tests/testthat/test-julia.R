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

skip_if_julia_unavailable <- function() {
  testthat::skip_if_not_installed("JuliaConnectoR")

  bindir <- find_julia_bindir()
  testthat::skip_if(bindir == "", "Julia executable not available.")

  Sys.setenv(JULIA_BINDIR = bindir)
  cs_ready <- tryCatch({
    JuliaConnectoR::juliaEval("using Circuitscape")
    TRUE
  }, error = function(e) FALSE)
  testthat::skip_if_not(cs_ready, "Circuitscape.jl is not available in Julia.")

  bindir
}

get_julia_pkg_data <- function(name) {
  data_env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "ResistanceGA2", envir = data_env)
  get(name, envir = data_env)
}

unwrap_julia_raster <- function(x) {
  if (inherits(x, "PackedSpatRaster")) {
    terra::unwrap(x)
  } else {
    x
  }
}

make_julia_example <- function(n = 5) {
  sample_pops <- get_julia_pkg_data("sample_pops")
  dc_list <- get_julia_pkg_data("Dc_list")
  pts <- terra::vect(sample_pops[[1]][seq_len(n), ], type = "points")
  gd <- ResistanceGA2::lower(dc_list[[1]][seq_len(n), seq_len(n)])
  keep <- rep(c(1L, 0L), length.out = length(gd))
  covariates <- data.frame(elev = seq_along(gd) / 100)

  list(
    pts = pts,
    gd = gd,
    keep = keep,
    covariates = covariates,
    response_cov = gd + covariates$elev
  )
}

test_that("Run_CS.jl auto-creates EXPORT.dir and returns a current map", {
  bindir <- skip_if_julia_unavailable()
  sample_pops <- get_julia_pkg_data("sample_pops")
  dc_list <- get_julia_pkg_data("Dc_list")
  raster_orig <- unwrap_julia_raster(get_julia_pkg_data("raster_orig"))

  pts <- terra::vect(sample_pops[[1]][1:5, ], type = "points")
  gd <- ResistanceGA2::lower(dc_list[[1]][1:5, 1:5])

  jl_inputs <- ResistanceGA2::jl.prep(
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
    ResistanceGA2::Run_CS.jl(
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

test_that("jl.prep handles pairs_to_include, covariates, and write.files", {
  bindir <- skip_if_julia_unavailable()
  example <- make_julia_example()

  write_dir <- tempfile("rga2-jl-write-")
  on.exit(unlink(write_dir, recursive = TRUE, force = TRUE), add = TRUE)

  jl_inputs <- suppressMessages(
    suppressWarnings(
      ResistanceGA2::jl.prep(
        n.Pops = 5,
        response = example$response_cov,
        CS_Point.File = example$pts,
        covariates = example$covariates,
        formula = response ~ elev,
        JULIA_HOME = bindir,
        pairs_to_include = example$keep,
        write.files = write_dir,
        run_test = FALSE,
        silent = TRUE
      )
    )
  )

  expect_true(dir.exists(write_dir))
  expect_true(file.exists(jl_inputs$CS_Point.File))
  expect_true(file.exists(jl_inputs$pairs_to_include))
  expect_identical(jl_inputs$keep, example$keep)
  expect_equal(nrow(jl_inputs$df), sum(example$keep))
  expect_equal(nrow(jl_inputs$covariates), sum(example$keep))
  expect_equal(ncol(jl_inputs$ZZ), sum(example$keep))
  expect_true(all(c("gd", "elev", "pop") %in% names(jl_inputs$df)))
  expect_match(paste(deparse(jl_inputs$formula), collapse = " "), "cd")
})

test_that("Run_CS.jl returns a full matrix with excluded pairs marked -1", {
  bindir <- skip_if_julia_unavailable()
  example <- make_julia_example()
  raster_orig <- unwrap_julia_raster(get_julia_pkg_data("raster_orig"))

  jl_inputs <- suppressWarnings(
    ResistanceGA2::jl.prep(
      n.Pops = 5,
      response = example$gd,
      CS_Point.File = example$pts,
      JULIA_HOME = bindir,
      pairs_to_include = example$keep,
      run_test = FALSE,
      silent = TRUE
    )
  )

  expect_warning(
    out <- ResistanceGA2::Run_CS.jl(
      jl.inputs = jl_inputs,
      r = raster_orig[["cont_orig"]],
      full.mat = TRUE
    ),
    "excluded pairs"
  )

  expect_equal(dim(out), c(5L, 5L))
  expect_true(any(out == -1, na.rm = TRUE))
  expect_true(all(diag(out) == 0))
})

test_that("Run_CS.ini writes INI overrides and returns full output", {
  bindir <- skip_if_julia_unavailable()
  example <- make_julia_example()
  raster_orig <- unwrap_julia_raster(get_julia_pkg_data("raster_orig"))

  export_dir <- tempfile("rga2-jl-ini-")
  on.exit(unlink(export_dir, recursive = TRUE, force = TRUE), add = TRUE)

  out <- suppressMessages(
    suppressWarnings(
      ResistanceGA2::Run_CS.ini(
        r = raster_orig[["cont_orig"]],
        CS_Point.File = example$pts,
        JULIA_HOME = bindir,
        EXPORT.dir = export_dir,
        return = "all",
        quiet = TRUE,
        connect_four_neighbors_only = TRUE,
        write_cur_maps = TRUE,
        write_cum_cur_map_only = TRUE,
        log_level = "critical",
        screenprint_log = FALSE,
        rm.files = FALSE
      )
    )
  )

  expect_true(file.exists(out$ini_file))
  expect_false(is.null(out$result))
  expect_true(is.matrix(out$pairwise_matrix))
  expect_equal(dim(out$pairwise_matrix), c(5L, 5L))
  expect_true(inherits(out$current_map, "SpatRaster"))

  ini_lines <- readLines(out$ini_file, warn = FALSE)
  ini_text <- paste(ini_lines, collapse = "\n")
  expect_match(ini_text, "scenario = pairwise", fixed = TRUE)
  expect_match(ini_text, "connect_four_neighbors_only = True", fixed = TRUE)
  expect_match(ini_text, "write_cur_maps = True", fixed = TRUE)
  expect_match(ini_text, "write_cum_cur_map_only = True", fixed = TRUE)
})

test_that("Run_CS.ini supports advanced-mode raster files supplied as terra objects", {
  bindir <- skip_if_julia_unavailable()
  raster_orig <- unwrap_julia_raster(get_julia_pkg_data("raster_orig"))

  habitat <- raster_orig[["cont_orig"]]
  source_r <- terra::rast(habitat)
  ground_r <- terra::rast(habitat)
  mask_r <- terra::rast(habitat)

  source_vals <- rep(0, terra::ncell(source_r))
  ground_vals <- rep(0, terra::ncell(ground_r))
  source_vals[5] <- 1
  ground_vals[terra::ncell(ground_r) - 4] <- 1
  terra::values(source_r) <- source_vals
  terra::values(ground_r) <- ground_vals
  terra::values(mask_r) <- 1

  export_dir <- tempfile("rga2-jl-advanced-")
  on.exit(unlink(export_dir, recursive = TRUE, force = TRUE), add = TRUE)

  out <- suppressMessages(
    suppressWarnings(
      ResistanceGA2::Run_CS.ini(
        r = habitat,
        JULIA_HOME = bindir,
        EXPORT.dir = export_dir,
        return = "all",
        quiet = TRUE,
        scenario = "advanced",
        source_file = source_r,
        ground_file = ground_r,
        ground_file_is_resistances = FALSE,
        mask_file = mask_r,
        write_cur_maps = TRUE,
        write_cum_cur_map_only = TRUE,
        rm.files = FALSE
      )
    )
  )

  expect_true(file.exists(out$ini_file))
  expect_false(is.null(out$result))
  expect_true(length(out$generated_files) > 0L)

  ini_lines <- readLines(out$ini_file, warn = FALSE)
  ini_text <- paste(ini_lines, collapse = "\n")
  expect_match(ini_text, "scenario = advanced", fixed = TRUE)
  expect_match(ini_text, "use_mask = True", fixed = TRUE)
  mask_line <- ini_lines[grepl("^mask_file = ", ini_lines)]
  source_line <- ini_lines[grepl("^source_file = ", ini_lines)]
  ground_line <- ini_lines[grepl("^ground_file = ", ini_lines)]
  expect_false(grepl("(Browse", mask_line, fixed = TRUE))
  expect_false(grepl("(Browse", source_line, fixed = TRUE))
  expect_false(grepl("(Browse", ground_line, fixed = TRUE))
})

test_that("SS_optim covers the Julia optimization branch", {
  bindir <- skip_if_julia_unavailable()
  example <- make_julia_example()
  raster_orig <- unwrap_julia_raster(get_julia_pkg_data("raster_orig"))

  jl_inputs <- ResistanceGA2::jl.prep(
    n.Pops = 5,
    response = example$gd,
    CS_Point.File = example$pts,
    JULIA_HOME = bindir,
    run_test = FALSE,
    silent = TRUE
  )
  ga_inputs <- ResistanceGA2::GA.prep(
    raster = raster_orig[["cont_orig"]],
    Results.dir = tempfile(pattern = "ss-optim-jl-"),
    pop.size = 10,
    maxiter = 1,
    run = 1,
    seed = 1,
    monitor = FALSE,
    quiet = TRUE
  )

  ss_out <- suppressMessages(
    suppressWarnings(
      ResistanceGA2::SS_optim(
        jl.inputs = jl_inputs,
        GA.inputs = ga_inputs,
        dist_mod = FALSE,
        null_mod = FALSE,
        diagnostic_plots = FALSE
      )
    )
  )

  expect_s3_class(ss_out, "resga_ss_optim")
  expect_s3_class(ss_out$AICc, "data.frame")
  expect_true("cont_orig" %in% ss_out$AICc$Surface)
})

test_that("MS_optim covers the Julia multisurface branch", {
  bindir <- skip_if_julia_unavailable()
  example <- make_julia_example()
  raster_orig <- unwrap_julia_raster(get_julia_pkg_data("raster_orig"))

  jl_inputs <- ResistanceGA2::jl.prep(
    n.Pops = 5,
    response = example$gd,
    CS_Point.File = example$pts,
    JULIA_HOME = bindir,
    run_test = FALSE,
    silent = TRUE
  )
  ga_inputs <- ResistanceGA2::GA.prep(
    raster = raster_orig[[c("cat_orig", "cont_orig")]],
    Results.dir = tempfile(pattern = "ms-optim-jl-"),
    pop.size = 10,
    maxiter = 1,
    run = 1,
    seed = 1,
    monitor = FALSE,
    quiet = TRUE
  )

  ms_out <- suppressMessages(
    suppressWarnings(
      ResistanceGA2::MS_optim(
        jl.inputs = jl_inputs,
        GA.inputs = ga_inputs,
        diagnostic_plots = FALSE
      )
    )
  )

  expect_s3_class(ms_out, "resga_ms_optim")
  expect_s3_class(ms_out$AICc.tab, "data.frame")
  expect_true(all(c("AICc", "LL", "R2m") %in% names(ms_out$AICc.tab)))
  expect_true(nrow(ms_out$percent.contribution) >= 2L)
})
