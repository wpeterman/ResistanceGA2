make_internal_optimizer_fixture <- function() {
  cont <- terra::rast(
    nrows = 15,
    ncols = 15,
    xmin = 0,
    xmax = 10,
    ymin = 0,
    ymax = 10
  )
  set.seed(7)
  terra::values(cont) <- runif(terra::ncell(cont), 1, 10)
  names(cont) <- "continuous"

  breaks <- stats::quantile(
    terra::values(cont),
    probs = c(1 / 3, 2 / 3),
    na.rm = TRUE
  )
  cat_r <- terra::classify(
    cont,
    rbind(
      c(-Inf, breaks[1], 1),
      c(breaks[1], breaks[2], 2),
      c(breaks[2], Inf, 3)
    )
  )
  names(cat_r) <- "categorical"

  coords <- matrix(
    c(1, 1,
      2, 8,
      5, 3,
      8, 9,
      9, 2),
    ncol = 2,
    byrow = TRUE
  )
  pts <- terra::vect(coords, type = "points")

  list(
    continuous = cont,
    categorical = cat_r,
    points = pts,
    base_gdist = ResistanceGA2::gdist.prep(
      n.Pops = nrow(coords),
      samples = pts,
      method = "costDistance"
    )
  )
}

get_internal_pkg_data <- function(name) {
  data_env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "ResistanceGA2", envir = data_env)
  get(name, envir = data_env)
}

unwrap_internal_raster <- function(x) {
  if (inherits(x, "PackedSpatRaster")) {
    terra::unwrap(x)
  } else {
    x
  }
}

test_that(".onLoad unwraps packed SpatRaster lazydata objects", {
  ns <- asNamespace("ResistanceGA2")
  lazydata_env <- ns[[".__NAMESPACE__."]][["lazydata"]]

  skip_if(is.null(lazydata_env), "ResistanceGA2 lazydata environment is unavailable.")

  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  terra::values(r) <- 1:4
  assign("codex_test_packed", terra::wrap(r), envir = lazydata_env)
  on.exit(rm("codex_test_packed", envir = lazydata_env), add = TRUE)

  ResistanceGA2:::.onLoad(NULL, "ResistanceGA2")

  expect_s4_class(get("codex_test_packed", envir = lazydata_env), "SpatRaster")
})

test_that("matrix expansion helpers preserve expected pair structure", {
  pop_n <- c(2, 1, 2)
  pop_mat <- matrix(
    c(0, 1, 2,
      1, 0, 3,
      2, 3, 0),
    nrow = 3,
    byrow = TRUE
  )

  keep <- ResistanceGA2:::expand.keep(pop_n)
  expanded_matrix <- ResistanceGA2:::expand.mat(pop_mat, pop_n, format = "matrix")
  expanded_vector <- ResistanceGA2:::expand.mat(pop_mat, pop_n)
  expanded_from_vector <- ResistanceGA2:::expand.mat(ResistanceGA2::lower(pop_mat), pop_n)
  expanded_full_vector <- ResistanceGA2:::expand.mat_vec(pop_n, pop_mat)

  expect_length(keep, choose(sum(pop_n), 2))
  expect_equal(sum(keep), 8)
  expect_equal(dim(expanded_matrix), c(sum(pop_n), sum(pop_n)))
  expect_length(expanded_vector, sum(keep))
  expect_equal(expanded_vector, expanded_from_vector)
  expect_equal(expanded_full_vector[keep == 1], expanded_vector)
})

test_that("expand.mat_ replicates from/to population pairs by sample counts", {
  pop_n <- c(2, 1, 2)
  from_mat <- matrix(
    c(0, 0, 0,
      1, 0, 0,
      1, 2, 0),
    nrow = 3,
    byrow = TRUE
  )
  to_mat <- matrix(
    c(0, 0, 0,
      2, 0, 0,
      3, 3, 0),
    nrow = 3,
    byrow = TRUE
  )

  out <- ResistanceGA2:::expand.mat_(from_mat, to_mat, pop_n)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 10L)
  expect_equal(as.character(out$pop1)[1:3], c("1", "1", "1"))
  expect_equal(as.character(out$pop2)[1:3], c("2", "2", "2"))
})

test_that("Resistance.Opt_AICc returns finite objectives for continuous and categorical inputs", {
  resistance_surfaces <- unwrap_internal_raster(get_internal_pkg_data("resistance_surfaces"))
  samples_df <- get_internal_pkg_data("samples")
  pts <- terra::vect(samples_df[, 2:3], type = "points")
  response <- ResistanceGA2::lower(as.matrix(dist(samples_df[, 2:3])))
  gdist_inputs <- ResistanceGA2::gdist.prep(
    n.Pops = nrow(samples_df),
    response = response,
    samples = pts,
    method = "costDistance"
  )

  cont_raster <- terra::subset(resistance_surfaces, "continuous")
  ga_cont <- suppressMessages(
    ResistanceGA2::GA.prep(
      raster = cont_raster,
      Results.dir = tempfile(pattern = "opt-aicc-cont-"),
      monitor = FALSE,
      quiet = TRUE
    )
  )
  parm_cont <- c(3, 2.5, 50)

  obj_min <- ResistanceGA2:::Resistance.Opt_AICc(
    PARM = parm_cont,
    Resistance = cont_raster,
    gdist.inputs = gdist_inputs,
    GA.inputs = ga_cont,
    Min.Max = "min",
    iter = 1,
    quiet = TRUE
  )
  obj_max <- ResistanceGA2:::Resistance.Opt_AICc(
    PARM = parm_cont,
    Resistance = cont_raster,
    gdist.inputs = gdist_inputs,
    GA.inputs = ga_cont,
    Min.Max = "max",
    iter = 1,
    quiet = TRUE
  )

  cat_raster <- terra::subset(resistance_surfaces, "categorical")
  ga_cat <- suppressMessages(
    ResistanceGA2::GA.prep(
      raster = cat_raster,
      Results.dir = tempfile(pattern = "opt-aicc-cat-"),
      monitor = FALSE,
      quiet = TRUE
    )
  )
  parm_cat <- c(1, 5, 10)
  cat_map <- data.frame(
    id = terra::unique(cat_raster)[[1]],
    parm = parm_cat
  )

  obj_cat <- ResistanceGA2:::Resistance.Opt_AICc(
    PARM = parm_cat,
    Resistance = cat_raster,
    gdist.inputs = gdist_inputs,
    GA.inputs = ga_cat,
    Min.Max = "min",
    iter = 1,
    quiet = TRUE
  )

  expect_true(is.finite(obj_min))
  expect_equal(obj_max, -obj_min, tolerance = 1e-8)
  expect_true(is.finite(obj_cat))
})

test_that("scaling, equation, and raster helpers return stable results", {
  varying <- c(2, 4, 6)
  constant <- c(5, 5, 5)

  expect_equal(ResistanceGA2:::OPTIM.DIRECTION("max"), -1)
  expect_equal(ResistanceGA2:::OPTIM.DIRECTION("min"), 1)
  expect_identical(
    names(ResistanceGA2:::Cont.Param(c(2.5, 50))),
    c("shape_opt", "max")
  )
  expect_equal(ResistanceGA2:::SCALE.vector(varying, 0, 1), c(0, 0.5, 1))
  expect_equal(ResistanceGA2:::SCALE.vector(constant, 1, 10), c(1, 1, 1))
  expect_equal(ResistanceGA2:::SCALE.vector(constant, 3, 3), c(0, 0, 0))

  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  terra::values(r) <- c(2, 4, 6, 8)
  scaled_r <- ResistanceGA2:::SCALE(r, 0, 1)
  expect_equal(range(terra::values(scaled_r), na.rm = TRUE), c(0, 1))

  const_r <- r
  terra::values(const_r) <- 5
  expect_equal(unique(as.vector(terra::values(ResistanceGA2:::SCALE(const_r, 1, 10)))), 1)
  expect_equal(unique(as.vector(terra::values(ResistanceGA2:::SCALE(const_r, 3, 3)))), 0)

  expect_equal(ResistanceGA2:::eq.set(list("A")), 1:9)
  expect_equal(ResistanceGA2:::eq.set(list("M")), c(1, 3, 5, 7, 9))
  expect_equal(ResistanceGA2:::eq.set(list("R")), c(2, 4, 6, 8, 9))
  expect_equal(ResistanceGA2:::get.EQ(3), "Monomolecular")
  expect_equal(ResistanceGA2:::get.EQ("Monomolecular"), 3)

  names(r) <- "helper_raster"
  expect_equal(sort(ResistanceGA2:::rast_unique(terra::classify(r, cbind(c(-Inf, 4, 6), c(4, 6, Inf), 1:3)))), 1:3)
})

test_that("sampling and matrix-reading helpers return valid structures", {
  set.seed(123)
  cat_starts <- ResistanceGA2:::sv.cat(levels = 3, pop.size = 4, min = 2, max = 10)
  cont_starts <- ResistanceGA2:::sv.cont.nG(
    direction = "Increase",
    pop.size = 5,
    max = 100,
    eqs = 1:9
  )

  expect_equal(dim(cat_starts), c(4, 3))
  expect_true(all(cat_starts[, 1] == 1))
  expect_equal(dim(cont_starts), c(5, 3))
  expect_true(all(cont_starts[, 1] %in% 1:9))
  expect_length(ResistanceGA2:::Increase.starts(1), 4L)
  expect_length(ResistanceGA2:::Increase.starts.nG(1), 3L)

  mat_file <- tempfile(fileext = ".txt")
  writeLines(
    c(
      "id 1 2 3",
      "1 0 1 2",
      "2 1 0 3",
      "3 2 3 0"
    ),
    con = mat_file
  )
  on.exit(unlink(mat_file), add = TRUE)

  expect_equal(ResistanceGA2:::read.matrix(mat_file), c(1, 2, 3))
  expect_equal(
    unname(as.matrix(ResistanceGA2:::read.matrix2(mat_file))),
    matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), nrow = 3, byrow = TRUE)
  )
})

test_that("ZZ helpers support expanded IDs, correlation groups, and column selection", {
  simple_id <- ResistanceGA2::To.From.ID(4)
  drop_simple <- rep(c(1L, 0L), length.out = nrow(simple_id))
  zz_select <- ResistanceGA2:::ZZ.mat_select(simple_id, drop = drop_simple)

  sp_loc <- terra::vect(
    cbind(c(0, 1, 10), c(0, 1, 10)),
    type = "points"
  )
  expect_warning(
    expanded_id <- ResistanceGA2::To.From.ID(
      sampled_pops = 3,
      pop_n = c(2, 1, 2),
      spLoc = sp_loc,
      nb = 2
    ),
    "sub-graphs"
  )
  drop_expanded <- rep(c(1L, 0L, 1L), length.out = nrow(expanded_id))
  zz_expanded <- ResistanceGA2:::ZZ.mat(expanded_id)
  zz_dropped <- ResistanceGA2:::ZZ.mat(expanded_id, drop = drop_expanded)

  expect_s4_class(zz_select, "dgCMatrix")
  expect_equal(ncol(zz_select), sum(drop_simple))
  expect_s4_class(zz_expanded, "dgCMatrix")
  expect_s4_class(zz_dropped, "dgCMatrix")
  expect_equal(ncol(zz_dropped), sum(drop_expanded))
})

test_that("Circuitscape file helpers materialize object inputs into temporary files", {
  habitat <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4)
  terra::values(habitat) <- seq_len(terra::ncell(habitat))
  names(habitat) <- "habitat"

  mask <- terra::rast(habitat)
  terra::values(mask) <- rep(c(1, 0), length.out = terra::ncell(mask))
  names(mask) <- "mask"

  source_pts <- terra::vect(
    matrix(c(0.5, 0.5,
             3.5, 3.5), ncol = 2, byrow = TRUE),
    type = "points"
  )
  terra::values(source_pts) <- data.frame(strength = c(2, 5))

  ground_tbl <- data.frame(
    value = c(10, 20),
    x = c(1.5, 2.5),
    y = c(1.5, 2.5)
  )

  config <- ResistanceGA2:::.cs_default_ini()
  config[["Mask file"]][["mask_file"]] <- mask
  config[["Options for advanced mode"]][["source_file"]] <- source_pts
  config[["Options for advanced mode"]][["ground_file"]] <- ground_tbl
  config[["Options for pairwise and one-to-all and all-to-one modes"]][["included_pairs_file"]] <-
    list(mode = "exclude", pairs = data.frame(a = c(1, 2), b = c(2, 3)))
  config[["Options for reclassification of habitat data"]][["reclass_file"]] <-
    data.frame(from = c(1, 2), to = c(10, 20))

  out <- ResistanceGA2:::.cs_materialize_file_options(config, scratch = tempdir())
  on.exit(unlink(out$temp_files, force = TRUE), add = TRUE)

  mask_path <- out$config[["Mask file"]][["mask_file"]]
  source_path <- out$config[["Options for advanced mode"]][["source_file"]]
  ground_path <- out$config[["Options for advanced mode"]][["ground_file"]]
  pairs_path <- out$config[["Options for pairwise and one-to-all and all-to-one modes"]][["included_pairs_file"]]
  reclass_path <- out$config[["Options for reclassification of habitat data"]][["reclass_file"]]

  expect_true(all(file.exists(c(mask_path, source_path, ground_path, pairs_path, reclass_path))))
  expect_match(readLines(pairs_path, warn = FALSE)[1], "mode exclude", fixed = TRUE)
  expect_equal(ncol(utils::read.table(source_path)), 3L)
  expect_equal(ncol(utils::read.table(ground_path)), 3L)
  expect_equal(ncol(utils::read.table(reclass_path)), 2L)
})

test_that("Circuitscape network value-file helper writes two-column tables", {
  network_values <- data.frame(node = c(1, 2, 5), value = c(1, 0, 3))
  out <- ResistanceGA2:::.cs_prepare_value_map_file(
    network_values,
    arg = "source_file",
    scratch = tempdir(),
    prefix = "circuitscape_network_",
    data_type = "network"
  )
  on.exit(unlink(out$temp_files, force = TRUE), add = TRUE)

  expect_true(file.exists(out$path))
  expect_equal(ncol(utils::read.table(out$path)), 2L)
})
