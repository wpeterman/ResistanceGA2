args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

find_project_root <- function() {
  if (length(file_arg) > 0) {
    raw_path <- sub("^--file=", "", file_arg[[1]])
    if (!identical(raw_path, "-") && file.exists(raw_path)) {
      script_path <- normalizePath(
        raw_path,
        winslash = "/",
        mustWork = TRUE
      )
      return(dirname(dirname(script_path)))
    }
  }

  cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(cwd, "DESCRIPTION"))) {
    return(cwd)
  }

  stop("Run this script from the package root or with Rscript.")
}

resolve_demo_profile <- function(profile = "default") {
  profile <- tolower(profile)
  if (!profile %in% c("fast", "default", "full")) {
    stop("`profile` must be one of 'fast', 'default', or 'full'.")
  }

  switch(
    profile,
    fast = list(pop.size = 10, maxiter = 2, run = 1),
    default = list(pop.size = 25, maxiter = 15, run = 5),
    full = list(pop.size = 50, maxiter = 50, run = 15)
  )
}

normalize_parallel_arg <- function(cores = NULL) {
  if (is.null(cores)) {
    return(FALSE)
  }

  cores <- suppressWarnings(as.integer(cores))
  if (is.na(cores) || cores < 1L) {
    stop("`cores` must be NULL or a positive integer.")
  }

  if (cores == 1L) FALSE else cores
}

build_ibe_ibr_demo_outputs <- function(profile = "default",
                                       cores = NULL,
                                       output_dir = NULL,
                                       seed = 123,
                                       quiet = TRUE) {
  project_root <- find_project_root()
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)

  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required to build the IBE/IBR demo outputs.")
  }

  pkgload::load_all(".", quiet = TRUE)
  library(terra)

  ga_settings <- resolve_demo_profile(profile)
  parallel_arg <- normalize_parallel_arg(cores)

  if (is.null(output_dir)) {
    output_dir <- file.path(project_root, "inst", "extdata", "ibe-ibr-demo")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  single_dir <- file.path(output_dir, "single_surface")
  multi_dir <- file.path(output_dir, "multi_surface")
  dir.create(single_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(multi_dir, recursive = TRUE, showWarnings = FALSE)

  data("rga2_demo", package = "ResistanceGA2", envir = environment())
  data("rga2_demo_covariates", package = "ResistanceGA2", envir = environment())

  if (inherits(rga2_demo_covariates, "PackedSpatRaster")) {
    rga2_demo_covariates <- terra::unwrap(rga2_demo_covariates)
  }

  r_covariate <- terra::subset(rga2_demo_covariates, c("err27", "ffp", "cti", "hli"))

  site_ids <- rga2_demo$sites$demo_id
  pair_matrix <- t(utils::combn(site_ids, 2))
  pair_index <- data.frame(
    from = pair_matrix[, 1],
    to = pair_matrix[, 2],
    stringsAsFactors = FALSE
  )

  coords <- rga2_demo$sample_coords
  site_env <- rga2_demo$site_environment

  geo <- as.numeric(scale(lower(as.matrix(dist(coords)))))
  ibe_diff <- site_pairwise_diff(site_env[, c("cti", "hli")], scale = TRUE)
  env_dist <- site_pairwise_dist(site_env[, c("cti", "hli")], scale = TRUE)

  predictor_cor <- stats::cor(
    data.frame(
      geo = geo,
      cti_diff = ibe_diff$cti_diff,
      hli_diff = ibe_diff$hli_diff,
      env_dist = env_dist$env_dist
    )
  )

  pair_data <- mlpe_data(
    response = rga2_demo$genetic_vector,
    pairs = pair_index,
    geo = geo,
    cti_diff = ibe_diff$cti_diff,
    hli_diff = ibe_diff$hli_diff,
    env_dist = env_dist$env_dist
  )

  baseline_models <- list(
    null = mlpe(
      response ~ 1 + (1 | pair),
      data = pair_data,
      pairs = c("from", "to")
    ),
    ibd = mlpe(
      response ~ geo + (1 | pair),
      data = pair_data,
      pairs = c("from", "to")
    ),
    ibe_cti = mlpe(
      response ~ cti_diff + (1 | pair),
      data = pair_data,
      pairs = c("from", "to")
    ),
    ibe_hli = mlpe(
      response ~ hli_diff + (1 | pair),
      data = pair_data,
      pairs = c("from", "to")
    ),
    ibe_both = mlpe(
      response ~ cti_diff + hli_diff + (1 | pair),
      data = pair_data,
      pairs = c("from", "to")
    ),
    ibd_ibe = mlpe(
      response ~ geo + cti_diff + hli_diff + (1 | pair),
      data = pair_data,
      pairs = c("from", "to")
    ),
    ibd_envdist = mlpe(
      response ~ geo + env_dist + (1 | pair),
      data = pair_data,
      pairs = c("from", "to")
    )
  )

  baseline_table <- do.call(
    rbind,
    lapply(names(baseline_models), function(name) {
      mod <- baseline_models[[name]]
      data.frame(
        model = name,
        AIC = AIC(mod),
        logLik = as.numeric(logLik(mod)),
        sigma = sigma(mod),
        row.names = NULL
      )
    })
  )
  baseline_table <- baseline_table[order(baseline_table$AIC), ]
  baseline_table$deltaAIC <- baseline_table$AIC - min(baseline_table$AIC)
  baseline_table$weight <- exp(-0.5 * baseline_table$deltaAIC)
  baseline_table$weight <- baseline_table$weight / sum(baseline_table$weight)

  baseline_fixef <- do.call(
    rbind,
    lapply(names(baseline_models), function(name) {
      cf <- lme4::fixef(baseline_models[[name]])
      data.frame(model = name, term = names(cf), estimate = unname(cf), row.names = NULL)
    })
  )

  pts <- terra::vect(
    coords,
    type = "points",
    crs = terra::crs(r_covariate, proj = TRUE)
  )
  workflow_covariates <- data.frame(
    geo = geo,
    cti_diff = ibe_diff$cti_diff,
    hli_diff = ibe_diff$hli_diff
  )

  gdist_inputs <- gdist.prep(
    n.Pops = nrow(coords),
    response = rga2_demo$genetic_vector,
    samples = pts,
    covariates = workflow_covariates,
    formula = gd ~ geo + cti_diff + hli_diff,
    method = "costDistance"
  )

  ibr_layers <- c("err27", "ffp")

  GA_inputs_ss <- do.call(
    GA.prep,
    c(
      list(
        raster = r_covariate[[ibr_layers]],
        Results.dir = single_dir,
        method = "LL",
        seed = seed,
        quiet = quiet,
        parallel = parallel_arg,
        keepBest = FALSE,
        monitor = FALSE
      ),
      ga_settings
    )
  )

  GA_inputs_ms <- do.call(
    GA.prep,
    c(
      list(
        raster = r_covariate[[ibr_layers]],
        Results.dir = multi_dir,
        method = "LL",
        seed = seed,
        quiet = quiet,
        parallel = parallel_arg,
        keepBest = FALSE,
        monitor = FALSE
      ),
      ga_settings
    )
  )

  ss_out <- SS_optim(
    gdist.inputs = gdist_inputs,
    GA.inputs = GA_inputs_ss,
    dist_mod = TRUE,
    null_mod = TRUE,
    diagnostic_plots = FALSE
  )

  ms_out <- MS_optim(
    gdist.inputs = gdist_inputs,
    GA.inputs = GA_inputs_ms,
    diagnostic_plots = FALSE
  )

  summary_bundle <- list(
    settings = list(
      profile = tolower(profile),
      parallel = parallel_arg,
      output_dir = output_dir,
      seed = seed,
      ibr_layers = ibr_layers,
      ga_settings = ga_settings,
      ibe_covariates = c("cti_diff", "hli_diff"),
      geographic_covariate = "geo"
    ),
    predictor_correlation = predictor_cor,
    baseline_table = baseline_table,
    baseline_fixef = baseline_fixef,
    ss_aicc = ss_out$AICc,
    ms_aicc = ms_out$AICc.tab,
    ms_percent_contribution = ms_out$percent.contribution
  )

  saveRDS(pair_data, file.path(output_dir, "pair_data.rds"))
  saveRDS(workflow_covariates, file.path(output_dir, "workflow_covariates.rds"))
  saveRDS(baseline_models, file.path(output_dir, "baseline_models.rds"))
  saveRDS(ss_out, file.path(output_dir, "single_surface_optim.rds"))
  saveRDS(ms_out, file.path(output_dir, "multi_surface_optim.rds"))
  saveRDS(summary_bundle, file.path(output_dir, "analysis_summary.rds"))

  utils::write.csv(
    baseline_table,
    file.path(output_dir, "baseline_model_table.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    baseline_fixef,
    file.path(output_dir, "baseline_fixed_effects.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    ss_out$AICc,
    file.path(output_dir, "single_surface_aicc.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    ms_out$AICc.tab,
    file.path(output_dir, "multi_surface_aicc.csv"),
    row.names = FALSE
  )

  message("Saved IBE/IBR demo outputs to: ", normalizePath(output_dir, winslash = "/"))

  invisible(
    list(
      output_dir = output_dir,
      summary_bundle = summary_bundle,
      ss_out = ss_out,
      ms_out = ms_out,
      baseline_models = baseline_models
    )
  )
}

if (sys.nframe() == 0L) {
  profile <- tolower(Sys.getenv("RGA2_DEMO_PROFILE", unset = "default"))
  legacy_fast <- identical(Sys.getenv("RGA2_DEMO_FAST", ""), "1")

  if (legacy_fast && profile == "default") {
    profile <- "fast"
  }

  cores_raw <- Sys.getenv("RGA2_DEMO_CORES", unset = "")
  cores <- if (nzchar(cores_raw)) as.integer(cores_raw) else NULL

  output_dir <- Sys.getenv("RGA2_DEMO_OUTPUT_DIR", unset = "")
  if (!nzchar(output_dir)) {
    output_dir <- NULL
  }

  build_ibe_ibr_demo_outputs(
    profile = profile,
    cores = cores,
    output_dir = output_dir,
    seed = 123,
    quiet = TRUE
  )
}
