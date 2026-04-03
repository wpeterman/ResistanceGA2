# =============================================================================
# generate_optimization_results.R
#
# Generates precomputed optimization results for use in vignettes.
#
# This script is run ONCE by the package developer (not during vignette
# rendering). Results are saved as .rds files in inst/extdata/vignette-results/
# and used by vignettes/ResistanceGA.Rmd via system.file().
#
# Also regenerates the package data objects stored in data/ (vignette_pairwise_
# outputs.rda and vignette_current_map.rda) that are produced by the earlier
# build_vignette_outputs.R. Run this script to refresh both sets of outputs.
#
# Usage (from the project root):
#   Rscript data-raw/generate_optimization_results.R
#
# Requirements:
#   - Julia must be installed with Circuitscape.jl
#   - JULIA_BINDIR environment variable must be set, or Julia must be on PATH
#   - pkgload must be installed
# =============================================================================
Sys.setenv(JULIA_BINDIR = "C:/Users/peterman.73/AppData/Local/Programs/Julia-1.11.4/bin/")

# Determine project root whether run via Rscript or sourced from RStudio
is_source_root <- function(path) {
  file.exists(file.path(path, "DESCRIPTION")) &&
    file.exists(file.path(path, "data-raw", "generate_optimization_results.R")) &&
    dir.exists(file.path(path, "data"))
}

args     <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  # Running via: Rscript data-raw/generate_optimization_results.R
  script_path  <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
  project_root <- dirname(dirname(script_path))
} else {
  # Running via source() inside RStudio: prefer the active project checkout.
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

  if (!is_source_root(project_root)) {
    project_root <- normalizePath(
      system.file(package = "ResistanceGA2"),
      winslash = "/", mustWork = FALSE
    )
  }

  if (!nzchar(project_root) || !is_source_root(project_root)) {
    stop(
      "Could not locate the ResistanceGA2 source tree. ",
      "Open the package project in RStudio, set the working directory to the ",
      "repository root, or run this script via Rscript from the checkout."
    )
  }
}

setwd(project_root)
message("Project root: ", project_root)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The 'pkgload' package is required. Install it with: install.packages('pkgload')")
}

# ---------------------------------------------------------------------------
# Load the package
# ---------------------------------------------------------------------------
message("\n[1/10] Loading ResistanceGA2 from source...")
pkgload::load_all(".", quiet = TRUE)
message("      Package loaded.")

required_data <- c("samples", "resistance_surfaces", "sample_pops", "Dc_list")
for (dataset_name in required_data) {
  data_file <- file.path(project_root, "data", paste0(dataset_name, ".rda"))
  if (!file.exists(data_file)) {
    stop("Required package data file not found: ", data_file)
  }
  load(data_file, envir = environment())
}

if (inherits(resistance_surfaces, "PackedSpatRaster")) {
  resistance_surfaces <- terra::unwrap(resistance_surfaces)
}

message("      Loaded package datasets: ", paste(required_data, collapse = ", "))

run_profile <- Sys.getenv("RGA2_OPTIM_PROFILE", unset = "full")
if (identical(run_profile, "smoke")) {
  ga_pop_size <- 5
  ga_maxiter <- 1
  ga_run <- 1
  ga_parallel <- FALSE
  all_comb_iters <- 5
} else if (identical(run_profile, "full")) {
  ga_pop_size <- 25
  ga_maxiter <- 200
  ga_run <- 25
  ga_parallel <- 8
  all_comb_iters <- 500
} else {
  stop("Unsupported RGA2_OPTIM_PROFILE='", run_profile, "'. Use 'full' or 'smoke'.")
}

message(
  "      Optimization profile: ", run_profile,
  " (pop.size=", ga_pop_size,
  ", maxiter=", ga_maxiter,
  ", run=", ga_run,
  ", parallel=", ga_parallel,
  ", all_comb iters=", all_comb_iters,
  ")"
)

# ---------------------------------------------------------------------------
# Julia discovery helper (copied from build_vignette_outputs.R)
# ---------------------------------------------------------------------------
find_julia_bindir <- function() {
  bindir <- Sys.getenv("JULIA_BINDIR", unset = "")
  if (nzchar(bindir) && dir.exists(bindir)) {
    return(normalizePath(bindir, winslash = "/", mustWork = TRUE))
  }

  candidates <- character()
  julia <- Sys.which("julia")
  if (nzchar(julia)) candidates <- c(candidates, julia)

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

# ---------------------------------------------------------------------------
# Output directories
# ---------------------------------------------------------------------------
vignette_results_dir <- file.path(project_root, "inst", "extdata", "vignette-results")
dir.create(vignette_results_dir, recursive = TRUE, showWarnings = FALSE)
message("\n[2/10] Output directory: ", vignette_results_dir)

# ---------------------------------------------------------------------------
# Set up data: small 25-sample dataset
# ---------------------------------------------------------------------------
message("\n[3/10] Setting up input data...")
set.seed(42)

pts_small <- terra::vect(samples[, 2:3], type = "points")
surface   <- terra::subset(resistance_surfaces, "continuous")

# Larger dataset for optimization (matches existing vignette examples)
pts_large <- terra::vect(sample_pops[[1]], type = "points")
gd_large  <- lower(Dc_list[[1]])

message("      Small dataset: ", nrow(samples), " populations")
message("      Large dataset: ", nrow(sample_pops[[1]]), " populations")

# ---------------------------------------------------------------------------
# [Part A] Regenerate vignette_pairwise_outputs and vignette_current_map
# (same outputs as build_vignette_outputs.R)
# ---------------------------------------------------------------------------
message("\n[4/10] Building gdistance vignette outputs (small dataset)...")

gdist_inputs_small <- gdist.prep(
  n.Pops  = nrow(samples),
  samples = pts_small,
  method  = "costDistance"
)

gdistance_commute <- Run_gdistance(
  gdist.inputs = gdist_inputs_small,
  r            = surface
)
if (inherits(gdistance_commute, "dist")) {
  gdistance_commute <- as.matrix(gdistance_commute)
}
message("      gdistance matrix: ", nrow(gdistance_commute), "x", ncol(gdistance_commute))

message("\n[5/10] Building Circuitscape vignette outputs (small dataset)...")

bindir <- find_julia_bindir()
Sys.setenv(JULIA_BINDIR = bindir)
message("      Julia bin dir: ", bindir)

tryCatch(
  JuliaConnectoR::juliaEval("using Circuitscape"),
  error = function(e) {
    stop("Circuitscape.jl is not available in Julia: ", conditionMessage(e))
  }
)
message("      Circuitscape.jl loaded.")

jl_inputs_small <- jl.prep(
  n.Pops        = nrow(samples),
  CS_Point.File = pts_small,
  JULIA_HOME    = bindir,
  run_test      = FALSE,
  silent        = TRUE,
  rm.files      = FALSE
)

circuitscape_matrix <- Run_CS.jl(
  jl.inputs  = jl_inputs_small,
  r          = surface,
  CurrentMap = FALSE,
  output     = "matrix"
)

export_dir <- tempfile("rga2-vignette-current-")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(export_dir, recursive = TRUE, force = TRUE), add = TRUE)

vignette_current_map_raw <- Run_CS.jl(
  jl.inputs  = jl_inputs_small,
  r          = surface,
  CurrentMap = TRUE,
  output     = "raster",
  EXPORT.dir = export_dir,
  rm.files   = FALSE
)

vignette_pairwise_outputs <- list(
  gdistance_commute   = gdistance_commute,
  circuitscape_matrix = circuitscape_matrix
)

save(
  vignette_pairwise_outputs,
  file     = file.path("data", "vignette_pairwise_outputs.rda"),
  compress = "bzip2"
)

vignette_current_map <- terra::wrap(vignette_current_map_raw)

save(
  vignette_current_map,
  file     = file.path("data", "vignette_current_map.rda"),
  compress = "bzip2"
)

message("      Saved data/vignette_pairwise_outputs.rda")
message("      Saved data/vignette_current_map.rda")

# ---------------------------------------------------------------------------
# [Part B] Create a "true" response surface with known parameters,
# then add noise to use as the genetic distance vector for optimization.
# This gives us a ground truth to verify the GA recovers.
#
# We combine the categorical and continuous layers from resistance_surfaces
# using known parameters, then simulate a noisy genetic distance vector.
#
# True parameters:
#   Categorical layer: class 1 = 1, class 2 = 250, class 3 = 75
#   Continuous layer:  transformation = Monomolecular (code 3), shape = 2.0, max = 100
#
# PARM vector for Combine_Surfaces() with two layers (cat + cont):
#   [cat_class1, cat_class2, cat_class3, trans_code, shape, max]
#   = c(1, 250, 75, 3, 2.0, 100)
# ---------------------------------------------------------------------------
message("\n[6/10] Creating synthetic response surface and genetic distances...")

# Use the resistance_surfaces stack (categorical + continuous) for demonstration
raster_two <- terra::subset(resistance_surfaces, c("categorical", "continuous"))

# GA settings for the true-response demo
ga_inputs_demo <- GA.prep(
  raster      = raster_two,
  Results.dir = vignette_results_dir,
  method      = "LL",
  pop.size    = ga_pop_size,
  maxiter     = ga_maxiter,
  run         = ga_run,
  parallel    = ga_parallel,
  monitor     = FALSE,
  quiet       = TRUE
)

gdist_inputs_demo <- gdist.prep(
  n.Pops   = nrow(samples),
  samples  = pts_small,
  method   = "costDistance"
)

# Known true parameters:
#   For categorical layer (3 classes): c1=1, c2=250, c3=75
#   For continuous layer: transform=Monomolecular (code 3), shape=2.0, max=100
true_parm <- c(1, 250, 75, 3, 2.0, 100)

# Generate the "true" composite surface
true_surface <- Combine_Surfaces(
  PARM         = true_parm,
  gdist.inputs = gdist_inputs_demo,
  GA.inputs    = ga_inputs_demo,
  rescale      = TRUE
)
message("      True composite surface created.")

# Compute true pairwise resistance distances
gdist_inputs_true <- gdist.prep(
  n.Pops   = nrow(samples),
  samples  = pts_small,
  method   = "costDistance"
)

gd_true <- Run_gdistance(
  gdist.inputs = gdist_inputs_true,
  r            = true_surface
)
gd_true_vec <- lower(as.matrix(gd_true))

# Add Gaussian noise (30% of SD)
set.seed(42)
noise     <- rnorm(length(gd_true_vec), mean = 0, sd = sd(gd_true_vec) * 0.3)
gd_noisy  <- gd_true_vec + noise
message("      Synthetic genetic distances created (n=", length(gd_noisy), ")")
message("      True LL signal + 30% noise added.")

# Rebuild gdist.inputs with the noisy genetic distances as response
gdist_inputs_optim <- gdist.prep(
  n.Pops   = nrow(samples),
  response = gd_noisy,
  samples  = pts_small,
  method   = "costDistance"
)

# ---------------------------------------------------------------------------
# [7/10] Single-surface optimization
# ---------------------------------------------------------------------------
message("\n[7/10] Running SS_optim()...")
message("      pop.size=", ga_pop_size, ", maxiter=", ga_maxiter, ", run=", ga_run)
message("      This may take several minutes...")

ss_out <- SS_optim(
  gdist.inputs     = gdist_inputs_optim,
  GA.inputs        = ga_inputs_demo,
  dist_mod         = TRUE,
  null_mod         = TRUE,
  diagnostic_plots = FALSE
)

saveRDS(ss_out, file.path(vignette_results_dir, "ss_optim_results.rds"))
message("      Saved ss_optim_results.rds")
message("      Top model by AICc:")
print(ss_out$AICc[which.min(ss_out$AICc$AICc), ])

# ---------------------------------------------------------------------------
# [8/10] Multi-surface optimization
# ---------------------------------------------------------------------------
message("\n[8/10] Running MS_optim()...")

ms_out <- MS_optim(
  gdist.inputs     = gdist_inputs_optim,
  GA.inputs        = ga_inputs_demo,
  diagnostic_plots = FALSE
)

saveRDS(ms_out, file.path(vignette_results_dir, "ms_optim_results.rds"))
message("      Saved ms_optim_results.rds")
message("      MS percent contributions:")
print(ms_out$percent.contribution)

# ---------------------------------------------------------------------------
# [9/10] All-combinations analysis
# ---------------------------------------------------------------------------
message("\n[9/10] Running all_comb()...")

ga_inputs_allcomb <- GA.prep(
  raster      = raster_two,
  Results.dir = "all.comb",
  method      = "LL",
  pop.size    = ga_pop_size,
  maxiter     = ga_maxiter,
  run         = ga_run,
  parallel    = ga_parallel,
  monitor     = FALSE,
  quiet       = TRUE
)

all_comb_results_dir <- file.path(vignette_results_dir, "all_comb")
dir.create(all_comb_results_dir, recursive = TRUE, showWarnings = FALSE)

combo_out <- all_comb(
  gdist.inputs    = gdist_inputs_optim,
  GA.inputs       = ga_inputs_allcomb,
  results.dir     = all_comb_results_dir,
  max.combination = 2,
  iters           = all_comb_iters,
  replicate       = 1
)

saveRDS(combo_out, file.path(vignette_results_dir, "all_comb_results.rds"))
message("      Saved all_comb_results.rds")
message("      All-combinations summary:")
print(combo_out$summary.table)

# ---------------------------------------------------------------------------
# [10/10] Save bootstrap results
# all_comb() runs Resist.boot() internally. The results are already in
# combo_out$boot.results. We save them separately for easy access in
# the ResistanceGA vignette.
# ---------------------------------------------------------------------------
message("\n[10/10] Saving bootstrap results from all_comb output...")

boot_out <- combo_out$boot.results
saveRDS(boot_out, file.path(vignette_results_dir, "boot_results.rds"))
message("      Saved boot_results.rds (from combo_out$boot.results)")
if (!is.null(boot_out)) {
  message("      Bootstrap summary:")
  print(as.data.frame(boot_out))
} else {
  message("      (boot.results was NULL — all_comb may not have run bootstrap)")
}

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
message("\n=============================================================")
message("All vignette results generated successfully.")
message("Files saved to: ", vignette_results_dir)
message("  - ss_optim_results.rds")
message("  - ms_optim_results.rds")
message("  - all_comb_results.rds")
message("  - boot_results.rds")
message("Package data updated:")
message("  - data/vignette_pairwise_outputs.rda")
message("  - data/vignette_current_map.rda")
message("")
message("True parameters used to generate synthetic response:")
message("  Categorical layer: class1=1, class2=250, class3=75")
message("  Continuous layer:  Monomolecular, shape=2.0, max=100")
message("  Noise: 30% of SD added (seed = 42)")
message("")
message("NOTE: boot_results.rds is extracted from all_comb output.")
message("  Resist.boot() runs automatically inside all_comb().")
message("  To run standalone bootstrap after SS_optim/MS_optim,")
message("  see the ResistanceGA vignette for correct argument syntax.")
message("=============================================================\n")
