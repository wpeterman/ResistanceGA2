# Internal helpers for worker-safe raster payloads ---------------------------

.rga_parallel_workers <- function(parallel) {
  if (inherits(parallel, "cluster")) {
    return(TRUE)
  }

  if (is.numeric(parallel) && length(parallel) == 1L && !is.na(parallel)) {
    return(parallel > 1)
  }

  isTRUE(parallel)
}

.rga_wrap_raster_for_parallel <- function(raster, parallel) {
  if (.rga_parallel_workers(parallel) && inherits(raster, "SpatRaster")) {
    return(terra::wrap(raster))
  }

  raster
}

.rga_materialize_raster <- function(raster) {
  if (inherits(raster, "PackedSpatRaster")) {
    return(terra::unwrap(raster))
  }

  raster
}

.rga_prepare_parallel_inputs <- function(GA.inputs) {
  if (!.rga_parallel_workers(GA.inputs$parallel)) {
    return(GA.inputs)
  }

  worker.inputs <- GA.inputs

  if (!is.null(worker.inputs$Resistance.stack)) {
    worker.inputs$Resistance.stack <-
      .rga_wrap_raster_for_parallel(worker.inputs$Resistance.stack, TRUE)
  }

  if (!is.null(worker.inputs$raster)) {
    worker.inputs$raster <-
      .rga_wrap_raster_for_parallel(worker.inputs$raster, TRUE)
  }

  if (!is.null(worker.inputs$inputs$raster)) {
    worker.inputs$inputs$raster <-
      .rga_wrap_raster_for_parallel(worker.inputs$inputs$raster, TRUE)
  }

  if (!is.null(worker.inputs$gdistance.approx.Resistance.stack)) {
    worker.inputs$gdistance.approx.Resistance.stack <-
      .rga_wrap_raster_for_parallel(
        worker.inputs$gdistance.approx.Resistance.stack,
        TRUE
      )
  }

  worker.inputs
}
