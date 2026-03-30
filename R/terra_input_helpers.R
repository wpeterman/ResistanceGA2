.validate_spatraster <- function(x, arg = "x", nlyr = NULL) {
  if (!inherits(x, "SpatRaster")) {
    stop("`", arg, "` must be a terra::SpatRaster.")
  }

  if (!is.null(nlyr) && terra::nlyr(x) != nlyr) {
    stop("`", arg, "` must have exactly ", nlyr, " layer", if (nlyr == 1L) "" else "s", ".")
  }

  x
}

.validate_spatvector_points <- function(x, arg = "x") {
  if (!inherits(x, "SpatVector")) {
    stop("`", arg, "` must be a terra::SpatVector of point locations.")
  }

  geom_type <- unique(tolower(as.character(terra::geomtype(x))))
  if (length(geom_type) != 1L || geom_type != "points") {
    stop("`", arg, "` must be a terra::SpatVector of points.")
  }

  x
}

.point_coords <- function(x, arg = "x") {
  x <- .validate_spatvector_points(x, arg = arg)
  terra::crds(x)
}

.gdistance_raster <- function(x, arg = "r") {
  x <- .validate_spatraster(x, arg = arg, nlyr = 1L)
  raster::raster(x)
}
