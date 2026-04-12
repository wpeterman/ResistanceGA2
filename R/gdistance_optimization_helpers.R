# Internal helpers for gdistance optimization approximation ------------------

.rga_prepare_gdistance_optimization_inputs <- function(GA.inputs, gdist.inputs) {
  if (!.rga_use_upstream_gdistance_approx(gdist.inputs)) {
    return(GA.inputs)
  }

  approx <- .gdistance_commute_approx_settings(gdist.inputs)
  out <- GA.inputs

  out$gdistance.approx.Resistance.stack <-
    .rga_aggregate_gdistance_optimization_stack(GA.inputs, approx$factor)
  out$gdistance.approx.factor <- approx$factor
  out$gdistance.approx.scale <- approx$scale
  out
}

.rga_use_upstream_gdistance_approx <- function(gdist.inputs) {
  if (is.null(gdist.inputs) ||
      !identical(gdist.inputs$method, "commuteDistance")) {
    return(FALSE)
  }

  approx <- .gdistance_commute_approx_settings(gdist.inputs)
  identical(approx$method, "aggregate") && approx$factor > 1L
}

.rga_has_upstream_gdistance_approx <- function(GA.inputs) {
  !is.null(GA.inputs$gdistance.approx.Resistance.stack)
}

.rga_aggregate_gdistance_optimization_stack <- function(GA.inputs, factor) {
  r <- .rga_materialize_raster(GA.inputs$Resistance.stack)
  layer_types <- GA.inputs$surface.type %||% rep("cont", terra::nlyr(r))

  layers <- lapply(seq_len(terra::nlyr(r)), function(i) {
    fun <- if (identical(layer_types[[i]], "cat")) {
      .rga_modal_value
    } else {
      mean
    }

    terra::aggregate(r[[i]], fact = factor, fun = fun, na.rm = TRUE)
  })

  out <- do.call(c, layers)
  names(out) <- names(r)
  out
}

.rga_modal_value <- function(x, ...) {
  x <- x[!is.na(x)]
  if (!length(x)) {
    return(NA_real_)
  }

  values <- sort(unique(x))
  values[which.max(tabulate(match(x, values)))]
}

.rga_optimization_resistance <- function(Resistance, GA.inputs, iter = NULL) {
  if (!.rga_has_upstream_gdistance_approx(GA.inputs) || is.null(iter)) {
    return(Resistance)
  }

  stack <- .rga_materialize_raster(GA.inputs$gdistance.approx.Resistance.stack)
  stack[[iter]]
}

.rga_gdistance_optimization_inputs <- function(GA.inputs) {
  if (!.rga_has_upstream_gdistance_approx(GA.inputs)) {
    return(GA.inputs)
  }

  out <- GA.inputs
  stack <- .rga_materialize_raster(GA.inputs$gdistance.approx.Resistance.stack)
  out$Resistance.stack <- stack
  out$raster <- stack
  out
}

.rga_run_gdistance_optimization <- function(gdist.inputs,
                                            r,
                                            GA.inputs,
                                            return.error.value = FALSE) {
  if (!.rga_has_upstream_gdistance_approx(GA.inputs)) {
    return(Run_gdistance(
      gdist.inputs,
      r,
      return.error.value = return.error.value
    ))
  }

  cd <- Run_gdistance(
    gdist.inputs,
    r,
    return.error.value = return.error.value,
    commute.approx = "none"
  )

  if (!identical(cd, -99999) && isTRUE(GA.inputs$gdistance.approx.scale)) {
    cd <- cd * (GA.inputs$gdistance.approx.factor^2)
  }

  cd
}

.rga_run_gdistance_exact <- function(gdist.inputs, r, ...) {
  Run_gdistance(gdist.inputs, r, commute.approx = "none", ...)
}
