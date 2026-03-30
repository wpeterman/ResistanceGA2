.onLoad <- function(libname, pkgname) {
 # Auto-unwrap PackedSpatRaster objects created by terra::wrap()
 # so that data(...) returns ready-to-use SpatRaster objects.
  ns <- asNamespace(pkgname)
  ld <- ns[[".__NAMESPACE__."]][["lazydata"]]
  if (!is.null(ld)) {
    for (nm in ls(ld)) {
      obj <- get(nm, envir = ld)
      if (inherits(obj, "PackedSpatRaster")) {
        assign(nm, terra::unwrap(obj), envir = ld)
      }
    }
  }
}
