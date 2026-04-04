args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 0) {
  stop("Run this script with Rscript so the project root can be determined.")
}

script_path <- normalizePath(
  sub("^--file=", "", file_arg[[1]]),
  winslash = "/",
  mustWork = TRUE
)
project_root <- dirname(dirname(script_path))
setwd(project_root)

landgen_sha <- "d5b9161c6d324878bae2f7bf285611b3a4eda1e6"
genetit_sha <- "c162b456b26f6f07aece0f2b619c22a6769da1f6"

download_asset <- function(url, destfile, mode = "wb") {
  utils::download.file(url, destfile, mode = mode, quiet = TRUE)
  destfile
}

tmp_dir <- tempfile("rga2-demo-")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

ralu_site_path <- download_asset(
  sprintf(
    "https://raw.githubusercontent.com/jeffreyevans/GeNetIt/%s/data/ralu.site.rda",
    genetit_sha
  ),
  file.path(tmp_dir, "ralu.site.rda")
)

covariates_path <- download_asset(
  sprintf(
    "https://raw.githubusercontent.com/jeffreyevans/GeNetIt/%s/inst/extdata/covariates.tif",
    genetit_sha
  ),
  file.path(tmp_dir, "covariates.tif")
)

genetic_path <- download_asset(
  sprintf(
    "https://raw.githubusercontent.com/hhwagner1/LandGenCourse/%s/inst/extdata/ralu_dc.csv",
    landgen_sha
  ),
  file.path(tmp_dir, "ralu_dc.csv")
)

load(ralu_site_path)

rga2_demo_covariates <- terra::rast(covariates_path)
rga2_demo_genetic <- as.matrix(utils::read.csv(genetic_path, check.names = FALSE))
storage.mode(rga2_demo_genetic) <- "double"
rownames(rga2_demo_genetic) <- colnames(rga2_demo_genetic)

selected_sites <- c(1:3, 5:6, 8, 11:12, 21, 23, 25:26)
site_id_map <- c(
  AirplaneLake = "Airplane",
  BachelorMeadow = "Bachelor",
  BarkingFoxLake = "BarkingFox",
  BobLake = "Bob",
  CacheLake = "Cache",
  EggWhiteLake = "Egg",
  FrogPondLake = "Frog",
  GentianLake = "GentianL",
  ParagonLake = "ParagonL",
  PotholeLake = "Pothole",
  ShipIslandLake = "ShipIsland",
  SkyhighLake = "Skyhigh"
)

site_vector <- terra::vect(ralu.site)[selected_sites]
site_data <- as.data.frame(site_vector, geom = "XY")
site_data$demo_id <- unname(site_id_map[site_data$SiteName])

if (anyNA(site_data$demo_id)) {
  stop("Failed to map one or more selected sites onto demo IDs.")
}

site_data <- site_data[, c("demo_id", setdiff(names(site_data), "demo_id"))]
rownames(site_data) <- site_data$demo_id

if (!identical(site_data$demo_id, colnames(rga2_demo_genetic))) {
  stop("Selected site order does not match the imported genetic matrix.")
}

site_environment <- terra::extract(rga2_demo_covariates, site_vector)
site_environment$ID <- NULL
site_environment <- data.frame(
  demo_id = site_data$demo_id,
  site_environment,
  row.names = site_data$demo_id,
  check.names = FALSE
)

sample_coords <- as.matrix(site_data[, c("x", "y")])
rownames(sample_coords) <- site_data$demo_id
colnames(sample_coords) <- c("x", "y")

rga2_demo <- list(
  sites = site_data,
  sample_coords = sample_coords,
  site_environment = site_environment,
  genetic_matrix = rga2_demo_genetic,
  genetic_vector = rga2_demo_genetic[lower.tri(rga2_demo_genetic)]
)

rga2_demo_covariates <- terra::wrap(rga2_demo_covariates)

save(
  rga2_demo,
  file = file.path("data", "rga2_demo.rda"),
  compress = "bzip2"
)

save(
  rga2_demo_covariates,
  file = file.path("data", "rga2_demo_covariates.rda"),
  compress = "bzip2"
)

save(
  rga2_demo_genetic,
  file = file.path("data", "rga2_demo_genetic.rda"),
  compress = "bzip2"
)

message(
  "Saved data/rga2_demo.rda, data/rga2_demo_covariates.rda, and ",
  "data/rga2_demo_genetic.rda"
)
