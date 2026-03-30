
get.name <- function(x) {
  nm <- deparse(substitute(x))
  return(nm)
}


Monomolecular <- function(r, parm) {
  parm[3] * (1 - exp(-1 * r / parm[2])) + 1   # Monomolecular
}

Inv.Monomolecular <- function(r, parm) {
  if (inherits(r, "SpatRaster")) {
    R <- parm[3] * (exp(-1 * r / parm[2]))
    R <- (R - terra::global(R, "min", na.rm = TRUE)[[1]]) + 1
  } else {
    R <- parm[3] * (exp(-1 * r / parm[2]))
    R <- (R - min(R, na.rm = TRUE)) + 1
  }
  R
}

Inv.Rev.Monomolecular <- function(r, parm) {
  if (inherits(r, "SpatRaster")) {
    rev.rast <- SCALE((-1 * r), 0, 10)
    Inv.Monomolecular(rev.rast, parm)
  } else {
    rev.rast <- SCALE.vector((-1 * r), 0, 10)
    Inv.Monomolecular(rev.rast, parm)
  }
}

Rev.Monomolecular <- function(r, parm) {
  if (inherits(r, "SpatRaster")) {
    rev.rast <- SCALE((-1 * r), 0, 10)
    Monomolecular(rev.rast, parm)
  } else {
    rev.rast <- SCALE.vector((-1 * r), 0, 10)
    Monomolecular(rev.rast, parm)
  }
}


Ricker <- function(r, parm) {
  parm[3] * r * exp(-1 * r / parm[2]) + 1   # Ricker
}

Inv.Ricker <- function(r, parm) {
  if (inherits(r, "SpatRaster")) {
    R <- (-1 * parm[3]) * r * exp(-1 * r / parm[2]) - 1
    R <- SCALE(R,
               MIN = abs(terra::global(R, "max", na.rm = TRUE)[[1]]),
               MAX = abs(terra::global(R, "min", na.rm = TRUE)[[1]]))
  } else {
    R <- (-1 * parm[3]) * r * exp(-1 * r / parm[2]) - 1
    R <- SCALE.vector(R, MIN = abs(max(R, na.rm = TRUE)), MAX = abs(min(R, na.rm = TRUE)))
  }
  R
}

Inv.Rev.Ricker <- function(r, parm) {
  if (inherits(r, "SpatRaster")) {
    rev.rast <- SCALE((-1 * r), 0, 10)
    Inv.Ricker(rev.rast, parm)
  } else {
    rev.rast <- SCALE.vector((-1 * r), 0, 10)
    Inv.Ricker(rev.rast, parm)
  }
}

Rev.Ricker <- function(r, parm) {
  if (inherits(r, "SpatRaster")) {
    rev.rast <- SCALE((-1 * r), 0, 10)
    Ricker(rev.rast, parm)
  } else {
    rev.rast <- SCALE.vector((-1 * r), 0, 10)
    Ricker(rev.rast, parm)
  }
}
