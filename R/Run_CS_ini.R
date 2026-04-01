# Helpers for low-level Circuitscape INI execution -------------------------

.cs_default_ini <- function(is_resistance = TRUE) {
  list(
    "Options for advanced mode" = list(
      ground_file_is_resistances = TRUE,
      remove_src_or_gnd = "rmvsrc",
      ground_file = "(Browse for a raster mask file)",
      use_unit_currents = FALSE,
      source_file = "(Browse for a raster mask file)",
      use_direct_grounds = FALSE
    ),
    "Mask file" = list(
      mask_file = "(Browse for a raster mask file)",
      use_mask = FALSE
    ),
    "Calculation options" = list(
      low_memory_mode = FALSE,
      parallelize = FALSE,
      max_parallel = 0L,
      solver = "cg+amg",
      print_timings = FALSE,
      preemptive_memory_release = FALSE,
      print_rusages = FALSE
    ),
    "Short circuit regions (aka polygons)" = list(
      polygon_file = "(Browse for a short-circuit region file)",
      use_polygons = FALSE
    ),
    "Options for one-to-all and all-to-one modes" = list(
      use_variable_source_strengths = FALSE,
      variable_source_file = "(Browse for a short-circuit region file)"
    ),
    "Output options" = list(
      set_null_currents_to_nodata = TRUE,
      set_focal_node_currents_to_zero = FALSE,
      set_null_voltages_to_nodata = TRUE,
      compress_grids = FALSE,
      write_volt_maps = FALSE,
      write_cur_maps = 0L,
      output_file = NULL,
      write_cum_cur_map_only = FALSE,
      log_transform_maps = FALSE,
      write_max_cur_maps = FALSE
    ),
    "Version" = list(
      version = "4.0.5"
    ),
    "Options for reclassification of habitat data" = list(
      reclass_file = "(Browse for file with reclassification data)",
      use_reclass_table = FALSE
    ),
    "Logging Options" = list(
      log_level = "critical",
      log_file = "None",
      profiler_log_file = "None",
      screenprint_log = FALSE
    ),
    "Options for pairwise and one-to-all and all-to-one modes" = list(
      included_pairs_file = "(Browse for a file with pairs to include or exclude)",
      use_included_pairs = FALSE,
      point_file = "(Browse for file with locations of focal points or regions)"
    ),
    "Connection scheme for raster habitat data" = list(
      connect_using_avg_resistances = FALSE,
      connect_four_neighbors_only = FALSE
    ),
    "Habitat raster or graph" = list(
      habitat_map_is_resistances = is_resistance,
      habitat_file = NULL
    ),
    "Circuitscape mode" = list(
      data_type = "raster",
      scenario = "pairwise",
      precision = "None"
    )
  )
}

.cs_template_keys <- function(config) {
  unlist(lapply(config, names), use.names = FALSE)
}

.cs_flat_config <- function(config) {
  values <- unlist(config, recursive = FALSE, use.names = FALSE)
  names(values) <- .cs_template_keys(config)
  values
}

.cs_placeholder_value <- function(x) {
  is.character(x) &&
    length(x) == 1L &&
    startsWith(x, "(") &&
    grepl("Browse for", x, fixed = TRUE)
}

.cs_ini_value <- function(key, value) {
  if (length(value) != 1L) {
    value <- paste(value, collapse = ",")
  }

  if (identical(key, "precision") && is.logical(value)) {
    return(if (isTRUE(value)) "single" else "None")
  }

  if (is.logical(value)) {
    return(if (isTRUE(value)) "True" else "False")
  }

  if (inherits(value, c("integer", "numeric"))) {
    return(as.character(value))
  }

  as.character(value)
}

.cs_apply_overrides <- function(config, overrides) {
  if (length(overrides) == 0L) {
    return(config)
  }

  if (is.null(names(overrides)) || any(!nzchar(names(overrides)))) {
    stop("All Circuitscape option overrides must be named.")
  }

  option_names <- gsub("\\.", "_", names(overrides))
  known_keys <- .cs_template_keys(config)
  unknown <- setdiff(option_names, known_keys)
  if (length(unknown) > 0L) {
    stop(
      "Unknown Circuitscape INI option(s): ",
      paste(unknown, collapse = ", "),
      ". Supply option names exactly as they appear in the Circuitscape .ini file."
    )
  }

  names(overrides) <- option_names
  for (i in seq_along(config)) {
    section_keys <- names(config[[i]])
    hits <- intersect(section_keys, option_names)
    if (length(hits) == 0L) {
      next
    }

    for (key in hits) {
      config[[i]][[key]] <- overrides[[key]]
    }
  }

  config
}

.cs_write_ini <- function(path, config) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  for (section in names(config)) {
    writeLines(paste0("[", section, "]"), con = con)
    for (key in names(config[[section]])) {
      value <- config[[section]][[key]]
      writeLines(
        paste0(key, " = ", .cs_ini_value(key, value)),
        con = con
      )
    }
    writeLines("", con = con)
  }

  invisible(path)
}

.cs_temp_dir <- function(scratch = NULL) {
  if (!is.null(scratch)) {
    return(paste0(normalizePath(scratch, winslash = "/", mustWork = TRUE), "/"))
  }

  if (Sys.info()[["sysname"]] == "Windows") {
    paste0(normalizePath(tempdir(), winslash = "/"), "/")
  } else {
    paste0(normalizePath(tempdir()), "/")
  }
}

.cs_normalize_output_dir <- function(EXPORT.dir = NULL, scratch = NULL) {
  if (is.null(EXPORT.dir)) {
    return(.cs_temp_dir(scratch = scratch))
  }

  EXPORT.dir <- paste0(
    sub("[/\\\\]+$", "", normalizePath(EXPORT.dir, winslash = "/", mustWork = FALSE)),
    "/"
  )
  if (!dir.exists(EXPORT.dir)) {
    dir.create(EXPORT.dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(EXPORT.dir)) {
      stop("Failed to create 'EXPORT.dir': ", EXPORT.dir)
    }
    message("Created 'EXPORT.dir': ", EXPORT.dir)
  }

  EXPORT.dir
}

.cs_prepare_habitat_file <- function(r,
                                     scratch = NULL,
                                     sanitize_raster = FALSE) {
  if (is.null(r)) {
    return(list(path = NULL, temp_files = character(), file_stub = "circuitscape"))
  }

  if (inherits(r, "SpatRaster")) {
    R <- .validate_spatraster(r, arg = "r", nlyr = 1L)
    if (isTRUE(sanitize_raster)) {
      mx.val <- terra::global(R, "max", na.rm = TRUE)[[1]]
      if (mx.val > 1e6) {
        R <- SCALE(R, 1, 1e6)
      }
      R <- terra::classify(R, matrix(c(-Inf, 0, 1), ncol = 3))
    } else {
      mn.val <- terra::global(R, "min", na.rm = TRUE)[[1]]
      if (is.finite(mn.val) && mn.val <= 0) {
        warning(
          "`r` contains values <= 0. Circuitscape raster inputs should be positive.",
          call. = FALSE
        )
      }
    }

    tmpdir <- .cs_temp_dir(scratch = scratch)
    temp_rast <- tempfile(pattern = "circuitscape_raster_", tmpdir = tmpdir, fileext = ".asc")
    terra::writeRaster(x = R, filename = temp_rast, overwrite = TRUE)

    stub <- names(R)
    if (!length(stub) || !nzchar(stub[[1]])) {
      stub <- "circuitscape"
    } else {
      stub <- stub[[1]]
    }

    return(list(
      path = normalizePath(temp_rast, winslash = "/", mustWork = TRUE),
      temp_files = normalizePath(temp_rast, winslash = "/", mustWork = TRUE),
      file_stub = stub
    ))
  }

  if (is.character(r) && length(r) == 1L) {
    if (!file.exists(r)) {
      stop("Raster file does not exist: ", r)
    }

    return(list(
      path = normalizePath(r, winslash = "/", mustWork = TRUE),
      temp_files = character(),
      file_stub = tools::file_path_sans_ext(basename(r))
    ))
  }

  stop("`r` must be a terra::SpatRaster, a path to an existing raster file, or NULL.")
}

.cs_prepare_point_file <- function(CS_Point.File, scratch = NULL) {
  if (is.null(CS_Point.File)) {
    return(list(path = NULL, temp_files = character(), n_points = NULL))
  }

  if (inherits(CS_Point.File, "SpatVector")) {
    pts <- .validate_spatvector_points(CS_Point.File, arg = "CS_Point.File")
    coords <- .point_coords(pts, arg = "CS_Point.File")
    site <- seq_len(nrow(coords))
    out <- data.frame(site, coords)

    tmpdir <- .cs_temp_dir(scratch = scratch)
    pt_file <- tempfile(pattern = "circuitscape_points_", tmpdir = tmpdir, fileext = ".txt")
    write.table(out, file = pt_file, col.names = FALSE, row.names = FALSE)

    return(list(
      path = normalizePath(pt_file, winslash = "/", mustWork = TRUE),
      temp_files = normalizePath(pt_file, winslash = "/", mustWork = TRUE),
      n_points = nrow(coords)
    ))
  }

  if (is.character(CS_Point.File) && length(CS_Point.File) == 1L) {
    if (!file.exists(CS_Point.File)) {
      stop("Point file does not exist: ", CS_Point.File)
    }

    n_points <- tryCatch(
      nrow(utils::read.table(CS_Point.File)),
      error = function(e) NULL
    )

    return(list(
      path = normalizePath(CS_Point.File, winslash = "/", mustWork = TRUE),
      temp_files = character(),
      n_points = n_points
    ))
  }

  stop("`CS_Point.File` must be a terra::SpatVector of points, a path to an existing point file, or NULL.")
}

.cs_temp_input_file <- function(prefix, scratch = NULL, fileext = ".txt") {
  normalizePath(
    tempfile(pattern = prefix, tmpdir = .cs_temp_dir(scratch), fileext = fileext),
    winslash = "/",
    mustWork = FALSE
  )
}

.cs_as_numeric_df <- function(x, arg) {
  if (is.matrix(x)) {
    x <- as.data.frame(x)
  }

  if (!is.data.frame(x)) {
    stop("`", arg, "` must be a data.frame, matrix, terra object, or file path.")
  }

  is_num <- vapply(x, function(col) is.numeric(col) || is.integer(col), logical(1))
  if (!all(is_num)) {
    stop("`", arg, "` must contain only numeric columns when supplied as a table.")
  }

  as.data.frame(x)
}

.cs_guess_value_col <- function(df, arg) {
  nm <- tolower(names(df))
  for (candidate in c("value", "strength", "current", "resistance",
                      "conductance", "ground", "source", "id")) {
    idx <- which(nm == candidate)
    if (length(idx) > 0L) {
      return(idx[[1]])
    }
  }

  numeric_cols <- which(vapply(df, is.numeric, logical(1)))
  if (length(numeric_cols) == 1L) {
    return(numeric_cols[[1]])
  }

  if (ncol(df) == 1L) {
    return(1L)
  }

  stop(
    "Could not determine which attribute column in `", arg,
    "` should be used as the Circuitscape value field."
  )
}

.cs_reorder_triples <- function(df, arg) {
  df <- .cs_as_numeric_df(df, arg)
  if (ncol(df) != 3L) {
    stop(
      "`", arg, "` must have exactly 3 columns in raster text-list form: value, x, y."
    )
  }

  nm <- tolower(names(df))
  value_idx <- match(TRUE, nm %in% c("value", "strength", "current",
                                     "resistance", "conductance", "ground",
                                     "source", "id"), nomatch = 0L)
  x_idx <- match(TRUE, nm %in% c("x", "lon", "longitude", "easting"), nomatch = 0L)
  y_idx <- match(TRUE, nm %in% c("y", "lat", "latitude", "northing"), nomatch = 0L)

  if (value_idx > 0L && x_idx > 0L && y_idx > 0L) {
    df <- df[, c(value_idx, x_idx, y_idx), drop = FALSE]
  }

  names(df) <- NULL
  df
}

.cs_write_table_file <- function(x,
                                 prefix,
                                 scratch = NULL,
                                 sep = "\t",
                                 col.names = FALSE,
                                 row.names = FALSE,
                                 quote = FALSE) {
  path <- .cs_temp_input_file(prefix, scratch = scratch, fileext = ".txt")
  write.table(
    x,
    file = path,
    sep = sep,
    col.names = col.names,
    row.names = row.names,
    quote = quote
  )

  list(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    temp_files = normalizePath(path, winslash = "/", mustWork = TRUE)
  )
}

.cs_write_raster_option <- function(x, arg, scratch = NULL, prefix = "circuitscape_opt_") {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      stop("File supplied to `", arg, "` does not exist: ", x)
    }

    return(list(
      path = normalizePath(x, winslash = "/", mustWork = TRUE),
      temp_files = character()
    ))
  }

  if (!inherits(x, "SpatRaster")) {
    stop("`", arg, "` must be a terra::SpatRaster or an existing file path.")
  }

  R <- .validate_spatraster(x, arg = arg, nlyr = 1L)
  path <- .cs_temp_input_file(prefix, scratch = scratch, fileext = ".asc")
  terra::writeRaster(x = R, filename = path, overwrite = TRUE)

  list(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    temp_files = normalizePath(path, winslash = "/", mustWork = TRUE)
  )
}

.cs_prepare_value_map_file <- function(x,
                                       arg,
                                       scratch = NULL,
                                       prefix = "circuitscape_values_",
                                       data_type = "raster") {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      stop("File supplied to `", arg, "` does not exist: ", x)
    }

    return(list(
      path = normalizePath(x, winslash = "/", mustWork = TRUE),
      temp_files = character()
    ))
  }

  data_type <- tolower(data_type)

  if (inherits(x, "SpatRaster")) {
    if (!identical(data_type, "raster")) {
      stop("`", arg, "` can only be supplied as a SpatRaster when `data_type = 'raster'`.")
    }
    return(.cs_write_raster_option(x, arg = arg, scratch = scratch, prefix = prefix))
  }

  if (inherits(x, "SpatVector")) {
    if (!identical(data_type, "raster")) {
      stop("`", arg, "` can only be supplied as a SpatVector when `data_type = 'raster'`.")
    }

    pts <- .validate_spatvector_points(x, arg = arg)
    attrs <- terra::values(pts)
    if (is.null(attrs) || ncol(attrs) == 0L) {
      stop("`", arg, "` SpatVector inputs must have at least one attribute column containing values.")
    }
    value_idx <- .cs_guess_value_col(attrs, arg)
    coords <- .point_coords(pts, arg = arg)
    table_out <- data.frame(
      value = attrs[[value_idx]],
      x = coords[, 1],
      y = coords[, 2]
    )

    return(.cs_write_table_file(
      table_out,
      prefix = prefix,
      scratch = scratch,
      sep = "\t",
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    ))
  }

  df <- .cs_as_numeric_df(x, arg)
  if (identical(data_type, "network")) {
    if (ncol(df) != 2L) {
      stop("`", arg, "` must have exactly 2 columns for Circuitscape network inputs: node_id, value.")
    }
    names(df) <- NULL
  } else {
    df <- .cs_reorder_triples(df, arg)
  }

  .cs_write_table_file(
    df,
    prefix = prefix,
    scratch = scratch,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
}

.cs_prepare_reclass_file <- function(x,
                                     arg = "reclass_file",
                                     scratch = NULL,
                                     prefix = "circuitscape_reclass_") {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      stop("File supplied to `", arg, "` does not exist: ", x)
    }

    return(list(
      path = normalizePath(x, winslash = "/", mustWork = TRUE),
      temp_files = character()
    ))
  }

  df <- .cs_as_numeric_df(x, arg)
  if (ncol(df) < 2L) {
    stop("`", arg, "` must have at least 2 columns when supplied as a table.")
  }
  names(df) <- NULL

  .cs_write_table_file(
    df,
    prefix = prefix,
    scratch = scratch,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
}

.cs_prepare_included_pairs_file <- function(x,
                                            arg = "included_pairs_file",
                                            scratch = NULL,
                                            prefix = "circuitscape_pairs_") {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      stop("File supplied to `", arg, "` does not exist: ", x)
    }

    return(list(
      path = normalizePath(x, winslash = "/", mustWork = TRUE),
      temp_files = character()
    ))
  }

  mode <- "include"
  pairs <- x

  if (is.list(x) && !is.data.frame(x) && !is.matrix(x)) {
    if (!all(c("mode", "pairs") %in% names(x))) {
      stop(
        "`", arg, "` lists must contain elements named `mode` and `pairs`."
      )
    }
    mode <- match.arg(tolower(as.character(x$mode)), c("include", "exclude"))
    pairs <- x$pairs
  } else if (is.data.frame(x) && ncol(x) == 2L) {
    nm <- tolower(names(x))
    if (identical(nm[[1]], "mode") && nm[[2]] %in% c("include", "exclude")) {
      mode <- nm[[2]]
    }
  }

  df <- .cs_as_numeric_df(pairs, arg)
  if (ncol(df) != 2L) {
    stop(
      "`", arg, "` must be a 2-column table of node pairs, a list with `mode` and `pairs`, or an existing file path."
    )
  }

  path <- .cs_temp_input_file(prefix, scratch = scratch, fileext = ".txt")
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste("mode", mode), con = con)
  write.table(df, file = con, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

  list(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    temp_files = normalizePath(path, winslash = "/", mustWork = TRUE)
  )
}

.cs_materialize_file_options <- function(config, scratch = NULL) {
  temp_files <- character()
  data_type <- tolower(as.character(config[["Circuitscape mode"]][["data_type"]]))

  value <- config[["Mask file"]][["mask_file"]]
  if (!is.null(value) && !.cs_placeholder_value(value) && !identical(value, "None")) {
    prep <- .cs_write_raster_option(value, arg = "mask_file", scratch = scratch, prefix = "circuitscape_mask_")
    config[["Mask file"]][["mask_file"]] <- prep$path
    temp_files <- c(temp_files, prep$temp_files)
  }

  value <- config[["Options for advanced mode"]][["source_file"]]
  if (!is.null(value) && !.cs_placeholder_value(value) && !identical(value, "None")) {
    prep <- .cs_prepare_value_map_file(
      value,
      arg = "source_file",
      scratch = scratch,
      prefix = "circuitscape_source_",
      data_type = data_type
    )
    config[["Options for advanced mode"]][["source_file"]] <- prep$path
    temp_files <- c(temp_files, prep$temp_files)
  }

  value <- config[["Options for advanced mode"]][["ground_file"]]
  if (!is.null(value) && !.cs_placeholder_value(value) && !identical(value, "None")) {
    prep <- .cs_prepare_value_map_file(
      value,
      arg = "ground_file",
      scratch = scratch,
      prefix = "circuitscape_ground_",
      data_type = data_type
    )
    config[["Options for advanced mode"]][["ground_file"]] <- prep$path
    temp_files <- c(temp_files, prep$temp_files)
  }

  value <- config[["Options for pairwise and one-to-all and all-to-one modes"]][["included_pairs_file"]]
  if (!is.null(value) && !.cs_placeholder_value(value) && !identical(value, "None")) {
    prep <- .cs_prepare_included_pairs_file(value, scratch = scratch)
    config[["Options for pairwise and one-to-all and all-to-one modes"]][["included_pairs_file"]] <- prep$path
    temp_files <- c(temp_files, prep$temp_files)
  }

  value <- config[["Options for reclassification of habitat data"]][["reclass_file"]]
  if (!is.null(value) && !.cs_placeholder_value(value) && !identical(value, "None")) {
    prep <- .cs_prepare_reclass_file(value, scratch = scratch)
    config[["Options for reclassification of habitat data"]][["reclass_file"]] <- prep$path
    temp_files <- c(temp_files, prep$temp_files)
  }

  list(config = config, temp_files = unique(temp_files))
}

.cs_output_base <- function(output_file) {
  output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
  sub("\\.out$", "", output_file)
}

.cs_existing_output_files <- function(output_base) {
  out_dir <- dirname(output_base)
  pattern <- paste0("^", basename(output_base))
  list.files(out_dir, pattern = pattern, full.names = TRUE, all.files = TRUE)
}

.cs_validate_file_options <- function(config) {
  keys_to_check <- c(
    "ground_file", "source_file", "mask_file", "polygon_file",
    "variable_source_file", "reclass_file", "included_pairs_file",
    "point_file", "habitat_file"
  )

  flat <- .cs_flat_config(config)
  for (key in intersect(keys_to_check, names(flat))) {
    value <- flat[[key]]
    if (is.null(value) || identical(value, "None") || .cs_placeholder_value(value)) {
      next
    }
    if (!file.exists(value)) {
      stop("Circuitscape option `", key, "` points to a file that does not exist: ", value)
    }
  }

  invisible(config)
}

.cs_validate_required_options <- function(config) {
  flat <- .cs_flat_config(config)
  scenario <- as.character(flat[["scenario"]])
  habitat_file <- flat[["habitat_file"]]

  if (is.null(habitat_file) || identical(habitat_file, "None") || .cs_placeholder_value(habitat_file)) {
    stop("A habitat raster / graph file is required. Supply `r` or override `habitat_file`.")
  }

  if (identical(flat[["data_type"]], "raster") &&
      scenario %in% c("pairwise", "one-to-all", "all-to-one") &&
      .cs_placeholder_value(flat[["point_file"]])) {
    stop(
      "A point file is required for Circuitscape raster scenario `",
      scenario,
      "`. Supply `CS_Point.File` or override `point_file`."
    )
  }

  if (identical(scenario, "advanced")) {
    if (.cs_placeholder_value(flat[["source_file"]])) {
      stop("For `scenario = 'advanced'`, supply `source_file`.")
    }
    if (.cs_placeholder_value(flat[["ground_file"]])) {
      stop("For `scenario = 'advanced'`, supply `ground_file`.")
    }
  }

  flag_requirements <- list(
    use_included_pairs = "included_pairs_file",
    use_mask = "mask_file",
    use_polygons = "polygon_file",
    use_reclass_table = "reclass_file",
    use_variable_source_strengths = "variable_source_file"
  )

  for (flag in names(flag_requirements)) {
    target <- flag_requirements[[flag]]
    if (isTRUE(flat[[flag]]) && .cs_placeholder_value(flat[[target]])) {
      stop(
        "Circuitscape option `", target, "` must be supplied when `",
        flag, " = TRUE`."
      )
    }
  }

  invisible(config)
}

.cs_extract_pairwise_matrix <- function(result, n_points = NULL) {
  mat <- tryCatch(as.matrix(result), error = function(e) NULL)
  if (is.null(mat) || length(dim(mat)) != 2L) {
    return(NULL)
  }

  out <- mat
  if (!is.null(n_points) && all(dim(mat) == c(n_points + 1L, n_points + 1L))) {
    out <- mat[-1, -1, drop = FALSE]
  }

  suppressWarnings(storage.mode(out) <- "double")
  out
}

.cs_load_current_map <- function(output_base) {
  cur_file <- paste0(output_base, "_cum_curmap.asc")
  if (!file.exists(cur_file)) {
    return(NULL)
  }

  terra::rast(cur_file)
}

.cs_point_count <- function(point_file) {
  if (is.null(point_file) || .cs_placeholder_value(point_file) || !file.exists(point_file)) {
    return(NULL)
  }

  tryCatch(
    nrow(utils::read.table(point_file)),
    error = function(e) NULL
  )
}

#' Run Circuitscape via Julia with direct INI-style options
#'
#' A low-level wrapper around \pkg{Circuitscape.jl} that writes a complete
#' Circuitscape \code{.ini} file from R, runs \code{compute()}, and returns the
#' raw result plus any generated output files.
#'
#' @param r A single-layer \code{SpatRaster} or a path to an existing raster
#'   file. If \code{NULL}, supply \code{habitat_file} via \code{cs_options} or
#'   \code{...}.
#' @param CS_Point.File A \code{\link[terra]{SpatVector}} of point locations or
#'   a path to an existing Circuitscape point file. Required for raster
#'   \code{pairwise}, \code{one-to-all}, and \code{all-to-one} scenarios unless
#'   \code{point_file} is overridden directly.
#' @param JULIA_HOME Path to the Julia \code{bin/} directory. If \code{NULL},
#'   \code{JuliaConnectoR} uses the current Julia configuration / environment.
#' @param EXPORT.dir Directory where the temporary \code{.ini}, \code{.out},
#'   and any Circuitscape output rasters are written.
#' @param mask_file Optional Circuitscape \code{mask_file}. May be a
#'   single-layer \code{SpatRaster} or an existing file path. When supplied,
#'   \code{use_mask} defaults to \code{TRUE} unless explicitly overridden.
#' @param source_file Optional Circuitscape \code{source_file}. May be a
#'   single-layer \code{SpatRaster}, a point \code{SpatVector} with a value
#'   attribute, a numeric table in Circuitscape text-list format, or an
#'   existing file path.
#' @param ground_file Optional Circuitscape \code{ground_file}. Same accepted
#'   forms as \code{source_file}. When using raster advanced mode and a raster
#'   is supplied, it is written as an ASCII grid.
#' @param included_pairs_file Optional Circuitscape
#'   \code{included_pairs_file}. May be an existing file path, a 2-column table
#'   of node pairs (defaults to \code{include} mode), or a list with elements
#'   \code{mode = "include"/"exclude"} and \code{pairs = <2-column table>}.
#'   When supplied, \code{use_included_pairs} defaults to \code{TRUE} unless
#'   explicitly overridden.
#' @param reclass_file Optional Circuitscape \code{reclass_file}. May be an
#'   existing file path or a numeric table, which is written to a temporary
#'   tab-delimited file as-is. When supplied, \code{use_reclass_table} defaults
#'   to \code{TRUE} unless explicitly overridden.
#' @param cs_options Named list of Circuitscape \code{.ini} option overrides.
#'   Names should match Circuitscape keys exactly (for example,
#'   \code{scenario}, \code{write_cur_maps}, \code{use_mask},
#'   \code{mask_file}).
#' @param return What to return: \code{"all"} (default) returns a list with the
#'   raw Julia result, parsed pairwise matrix when available, current map when
#'   available, the final config, and file paths; \code{"compute"} returns the
#'   raw Julia result; \code{"matrix"} returns a parsed pairwise matrix when
#'   available; \code{"current_map"} returns the cumulative current map raster.
#' @param rm.files Logical. Remove generated \code{.ini}, temporary inputs, and
#'   Circuitscape outputs after the run? Default = \code{FALSE}.
#' @param quiet Logical. Suppress Circuitscape progress output emitted through
#'   Julia? Default = \code{TRUE}.
#' @param scratch Optional scratch directory for temporary raster / point files.
#' @param is_resistance Logical. Should the input habitat map be treated as
#'   resistances (\code{TRUE}) or conductances (\code{FALSE}) by default?
#' @param sanitize_raster Logical. If \code{TRUE}, apply the same raster
#'   sanitizing used by \code{\link{Run_CS.jl}} before writing the temporary
#'   ASCII file. Default = \code{FALSE}.
#' @param ... Additional Circuitscape \code{.ini} option overrides supplied as
#'   named arguments. These are merged with \code{cs_options}.
#'
#' @return Depending on \code{return}, either the raw Julia result, a parsed
#'   pairwise matrix, a cumulative current-map raster, or a list containing all
#'   of those plus the written \code{.ini} file path, generated output files,
#'   and the final Circuitscape configuration used for the run.
#'
#' @details
#' This function complements \code{\link{Run_CS.jl}}. Whereas
#' \code{Run_CS.jl} is tailored to the package optimization workflow,
#' \code{Run_CS.ini} is intended for direct Circuitscape runs where the user
#' wants to specify Circuitscape settings in a way that closely mirrors the
#' \code{.ini} file.
#'
#' Option names in \code{cs_options} and \code{...} should use the same keys
#' that appear in a Circuitscape \code{.ini} file. The function currently
#' exposes the full raster-mode option template written by the package for
#' Circuitscape 4.0.5 and validates that overridden keys are known.
#'
#' The file-oriented arguments (\code{mask_file}, \code{source_file},
#' \code{ground_file}, \code{included_pairs_file}, and \code{reclass_file})
#' are convenience layers over the corresponding Circuitscape options. They are
#' especially useful when working directly with terra objects or small R tables.
#' For any more specialized file layout, users can still supply a file path
#' directly through these arguments or through \code{cs_options}.
#'
#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#'
#' @examples
#' \dontrun{
#' pts <- terra::vect(sample_pops[[1]], type = "points")
#'
#' cs.out <- Run_CS.ini(
#'   r = raster_orig[["cont_orig"]],
#'   CS_Point.File = pts,
#'   JULIA_HOME = Sys.getenv("JULIA_BINDIR"),
#'   scenario = "pairwise",
#'   connect_four_neighbors_only = TRUE,
#'   write_cur_maps = TRUE,
#'   write_cum_cur_map_only = TRUE,
#'   rm.files = FALSE
#' )
#' }
Run_CS.ini <- function(r = NULL,
                       CS_Point.File = NULL,
                       JULIA_HOME = NULL,
                       EXPORT.dir = NULL,
                       mask_file = NULL,
                       source_file = NULL,
                       ground_file = NULL,
                       included_pairs_file = NULL,
                       reclass_file = NULL,
                       cs_options = list(),
                       return = c("all", "compute", "matrix", "current_map"),
                       rm.files = FALSE,
                       quiet = TRUE,
                       scratch = NULL,
                       is_resistance = TRUE,
                       sanitize_raster = FALSE,
                       ...) {
  return <- match.arg(return)

  if (!is.list(cs_options)) {
    stop("`cs_options` must be a named list of Circuitscape option overrides.")
  }

  if (!is.null(JULIA_HOME)) {
    if (!dir.exists(JULIA_HOME)) {
      stop("Specified JULIA_HOME directory does not exist: ", JULIA_HOME)
    }
    Sys.setenv(JULIA_BINDIR = normalizePath(JULIA_HOME, winslash = "/", mustWork = TRUE))
  }
  JuliaConnectoR::juliaEval("using Circuitscape")

  export_dir <- .cs_normalize_output_dir(EXPORT.dir = EXPORT.dir, scratch = scratch)
  habitat <- .cs_prepare_habitat_file(
    r = r,
    scratch = scratch,
    sanitize_raster = sanitize_raster
  )
  points <- .cs_prepare_point_file(CS_Point.File = CS_Point.File, scratch = scratch)

  file_stub <- habitat$file_stub
  if (!nzchar(file_stub)) {
    file_stub <- "circuitscape"
  }

  ini_file <- normalizePath(
    file.path(export_dir, paste0(file_stub, ".ini")),
    winslash = "/",
    mustWork = FALSE
  )
  output_file <- normalizePath(
    file.path(export_dir, paste0(file_stub, ".out")),
    winslash = "/",
    mustWork = FALSE
  )

  config <- .cs_default_ini(is_resistance = is_resistance)
  if (!is.null(habitat$path)) {
    config[["Habitat raster or graph"]][["habitat_file"]] <- habitat$path
  }
  if (!is.null(points$path)) {
    config[["Options for pairwise and one-to-all and all-to-one modes"]][["point_file"]] <- points$path
  }
  config[["Output options"]][["output_file"]] <- output_file

  overrides <- c(
    cs_options[!vapply(cs_options, is.null, logical(1))],
    list(...))
  overrides <- overrides[!vapply(overrides, is.null, logical(1))]
  config <- .cs_apply_overrides(config, overrides)

  if (!is.null(mask_file)) {
    config[["Mask file"]][["mask_file"]] <- mask_file
    config[["Mask file"]][["use_mask"]] <- TRUE
  }
  if (!is.null(source_file)) {
    config[["Options for advanced mode"]][["source_file"]] <- source_file
  }
  if (!is.null(ground_file)) {
    config[["Options for advanced mode"]][["ground_file"]] <- ground_file
  }
  if (!is.null(included_pairs_file)) {
    config[["Options for pairwise and one-to-all and all-to-one modes"]][["included_pairs_file"]] <- included_pairs_file
    config[["Options for pairwise and one-to-all and all-to-one modes"]][["use_included_pairs"]] <- TRUE
  }
  if (!is.null(reclass_file)) {
    config[["Options for reclassification of habitat data"]][["reclass_file"]] <- reclass_file
    config[["Options for reclassification of habitat data"]][["use_reclass_table"]] <- TRUE
  }

  materialized <- .cs_materialize_file_options(config, scratch = scratch)
  config <- materialized$config

  .cs_validate_required_options(config)
  .cs_validate_file_options(config)

  output_file_final <- config[["Output options"]][["output_file"]]
  out_dir_final <- dirname(normalizePath(output_file_final, winslash = "/", mustWork = FALSE))
  if (!dir.exists(out_dir_final)) {
    dir.create(out_dir_final, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(out_dir_final)) {
      stop("Failed to create output directory for `output_file`: ", out_dir_final)
    }
  }

  .cs_write_ini(ini_file, config)

  raw_result <- NULL
  if (isTRUE(quiet)) {
    invisible(capture.output(
      raw_result <- JuliaConnectoR::juliaCall("compute", ini_file),
      type = "message"
    ))
  } else {
    raw_result <- JuliaConnectoR::juliaCall("compute", ini_file)
  }

  output_base <- .cs_output_base(config[["Output options"]][["output_file"]])
  point_count <- points$n_points
  if (is.null(point_count)) {
    point_count <- .cs_point_count(
      config[["Options for pairwise and one-to-all and all-to-one modes"]][["point_file"]]
    )
  }
  pairwise_matrix <- .cs_extract_pairwise_matrix(raw_result, n_points = point_count)
  current_map <- .cs_load_current_map(output_base)
  generated_files <- .cs_existing_output_files(output_base)
  temp_files <- unique(c(habitat$temp_files, points$temp_files, materialized$temp_files, ini_file))

  result <- switch(
    return,
    all = list(
      result = raw_result,
      pairwise_matrix = pairwise_matrix,
      current_map = current_map,
      ini_file = ini_file,
      output_dir = export_dir,
      output_base = output_base,
      generated_files = generated_files,
      config = config,
      input_files = list(
        habitat_file = config[["Habitat raster or graph"]][["habitat_file"]],
        point_file = config[["Options for pairwise and one-to-all and all-to-one modes"]][["point_file"]]
      )
    ),
    compute = raw_result,
    matrix = {
      if (is.null(pairwise_matrix)) {
        stop(
          "A pairwise matrix could not be parsed from the Circuitscape result. ",
          "Use `return = 'compute'` or `return = 'all'` for the raw Julia output."
        )
      }
      pairwise_matrix
    },
    current_map = {
      if (is.null(current_map)) {
        stop(
          "No cumulative current map was found. ",
          "Set options such as `write_cur_maps = TRUE` and ",
          "`write_cum_cur_map_only = TRUE` before requesting `return = 'current_map'`."
        )
      }
      current_map
    }
  )

  if (isTRUE(rm.files)) {
    cleanup_files <- unique(c(generated_files, temp_files))
    cleanup_files <- cleanup_files[file.exists(cleanup_files)]
    invisible(sapply(cleanup_files, unlink, force = TRUE))
  }

  result
}
