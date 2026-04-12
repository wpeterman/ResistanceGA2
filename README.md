# ResistanceGA2

**Optimize Resistance Surfaces Using Genetic Algorithms**

`ResistanceGA2` is a major revision of the original `ResistanceGA` package. It retains the core GA-based resistance surface optimization framework while modernizing the spatial infrastructure, expanding the analytical capabilities, and removing deprecated workflows. Both continuous and categorical resistance surfaces can be optimized individually (`SS_optim`) or jointly as a composite (`MS_optim`). All-combinations enumeration (`all_comb`) and bootstrap resampling (`Resist.boot`) are also supported.

---

## Installation

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("wpeterman/ResistanceGA2", build_vignettes = TRUE)

library(ResistanceGA2)
```

For Julia/Circuitscape support, see the `Julia_Guide` vignette after installation.

---

## What Is New in ResistanceGA2

The sections below describe every substantive change from the original `ResistanceGA`. Users migrating from `ResistanceGA` should read this section carefully before updating scripts.

### 1. terra replaces raster and sp

All spatial inputs and outputs now use the `terra` package. The legacy `raster` and `sp` classes are no longer accepted as direct inputs.

| Old (ResistanceGA)             | New (ResistanceGA2)                  |
|-------------------------------|--------------------------------------|
| `RasterLayer`, `RasterStack`  | `SpatRaster`                         |
| `SpatialPoints`, `SpatialPointsDataFrame` | `SpatVector` (point layer)  |
| `stack()`, `raster()`         | `terra::rast()`, `terra::vect()`     |
| `raster::subset()`, `dropLayer()` | `terra::subset()`, `[[i]]`       |

**Bundled datasets** have all been converted from serialized raster/sp S4 objects to wrapped `SpatRaster` objects (via `terra::wrap()`). An `.onLoad` hook automatically unwraps them when the package is attached, so `resistance_surfaces` and related objects are immediately usable as `SpatRaster` without any manual step.

```r
# Old ResistanceGA
library(raster)
r <- stack("surface.tif")

# ResistanceGA2
library(terra)
r <- rast("surface.tif")
```

### 2. JuliaConnectoR replaces JuliaCall

The Julia/Circuitscape interface now uses `JuliaConnectoR` instead of `JuliaCall`. The user-facing function signatures are unchanged (`jl.prep()`, `Run_CS.jl()`), but the underlying R-to-Julia bridge is different. If you have old scripts that call `JuliaCall` functions directly, update them to `JuliaConnectoR` equivalents.

Julia console output (the verbose `@info` and `@warn` messages from Circuitscape.jl) is now automatically suppressed during optimization runs. Set `silent = TRUE` in `jl.prep()` to suppress it during the initial test run as well.

### 3. Standalone Circuitscape .exe workflow removed

`CS.prep()` and `Run_CS()` (the functions that called the Circuitscape standalone executable) have been removed. Circuitscape runs are handled exclusively through Julia via `jl.prep()` and `Run_CS.jl()`.

| Old (ResistanceGA)  | New (ResistanceGA2)  |
|--------------------|----------------------|
| `CS.prep()`        | `jl.prep()`          |
| `Run_CS()`         | `Run_CS.jl()`        |
| `SS_optim(CS.inputs = ...)` | `SS_optim(jl.inputs = ...)` |

### 4. Grid search removed

`Grid.Search()` has been removed. The `akima` package dependency has been dropped. If you need a grid search over discrete parameter values, use the GA with a small population and restricted transformation families instead.

### 5. k.smooth removed

The `k.smooth` parameter used in some legacy transformation workflows has been removed. Transformations are now specified entirely through `select.trans` in `GA.prep()`.

### 6. Covariates in optimization (IBE-conditional IBR)

`gdist.prep()` and `jl.prep()` now accept `covariates` and `formula` arguments. When these are supplied, every candidate resistance surface during the GA search is evaluated in an MLPE model that already conditions on the specified fixed effects (e.g., geographic distance, environmental mismatch). This makes it possible to test whether a landscape resistance hypothesis explains genetic differentiation *after* accounting for IBD or IBE.

```r
# Condition resistance optimization on geographic distance and an
# environmental mismatch term
gdist_inputs <- gdist.prep(
  n.Pops     = nrow(coords),
  response   = gd_vector,
  samples    = pts,
  covariates = data.frame(geo = geo_dist, cti_diff = env_diff$cti_diff),
  formula    = gd ~ geo + cti_diff,
  method     = "costDistance"
)
```

See the `IBE_IBR_Framework` vignette for a full worked example.

### 7. New MLPE interface: mlpe() and mlpe_data()

Two new functions provide a cleaner, dyad-aware interface for fitting MLPE models directly (outside of the optimization workflow):

- **`mlpe_data(response, ..., pairs, labels)`**: Converts pairwise responses and named covariate vectors into a long-form data frame ready for `mlpe()`. Inputs can be lower-triangle vectors, symmetric matrices, or `dist` objects.
- **`mlpe(formula, data, pairs, ...)`**: Fits an MLPE mixed model using `lme4` syntax. The dyadic random-effect covariance structure is built automatically from the endpoint columns specified in `pairs`.

```r
pair_data <- mlpe_data(
  response = gd_vector,
  pairs    = pair_index,      # data frame with 'from' / 'to' columns
  geo      = geo_distance,
  cti_diff = env_mismatch$cti_diff
)

m_ibr <- mlpe(
  response ~ geo + cti_diff + cd + (1 | pair),
  data  = pair_data,
  pairs = c("from", "to")
)
```

The older `MLPE.lmm()`, `mlpe_rga()`, and `To.From.ID()` functions remain available for backward compatibility and are still used internally during optimization, but `mlpe()` / `mlpe_data()` is the recommended API for new standalone analyses.

### 8. IBE helper functions

Two new functions construct pairwise environmental covariates from site-level measurements, aligned to the lower-triangle ordering used throughout the package:

- **`site_pairwise_diff(data, scale, transform, ...)`**: Computes absolute (or squared) pairwise differences in site-level environmental variables — the standard representation for isolation-by-environment (IBE) covariates.
- **`site_pairwise_dist(data, scale, ...)`**: Computes multivariate Euclidean distance among sites in environmental space — an alternative single-term IBE representation.

```r
ibe_diff <- site_pairwise_diff(site_env[, c("cti", "hli")], scale = TRUE)
env_dist  <- site_pairwise_dist(site_env[, c("cti", "hli")], scale = TRUE)
```

### 9. S3 classes for optimization results

Optimization functions (`SS_optim`, `MS_optim`, `all_comb`) now return objects with S3 classes and associated `print()`, `summary()`, and `plot()` methods. The underlying list structure (with `$AICc`, `$cd`, `$MLPE.list`, etc.) is unchanged, so existing code continues to work. The methods add a cleaner interactive layer.

```r
ss_out <- SS_optim(gdist.inputs = gdist_inputs, GA.inputs = ga_inputs)

class(ss_out)
print(ss_out)
summary(ss_out)
plot(ss_out)
```

### 10. Low-level Circuitscape INI wrapper

A new `Run_CS.ini()` function provides direct control over Circuitscape run parameters through the `.ini` file interface. This is useful when you need features not exposed by `Run_CS.jl()`, such as advanced mode (source/ground inputs), mask files, short-circuit regions, or included-pairs files.

```r
cs_out <- Run_CS.ini(
  r                        = surface,
  CS_Point.File            = pts,
  JULIA_HOME               = JULIA_HOME,
  scenario                 = "pairwise",
  connect_four_neighbors_only = TRUE,
  write_cur_maps           = TRUE,
  write_cum_cur_map_only   = TRUE,
  return                   = "all",
  rm.files                 = FALSE
)
```

See the `Circuitscape_INI` vignette for full documentation.

### 11. GA.prep Results.dir behavior

`GA.prep()` now actively manages the `Results.dir` path:

- If the specified directory does not exist and R is running interactively, the user is prompted to create it.
- In non-interactive sessions, the directory is created automatically.
- If the directory already exists and is not empty, a message warns that existing results may be overwritten.

### 12. Idaho demo dataset

A new Idaho landscape genetics demo dataset ships with the package:

- **`rga2_demo`**: A list containing site coordinates (`sample_coords`), site-level environmental variables (`site_environment`), a genetic distance vector (`genetic_vector`), and site metadata (`sites`).
- **`rga2_demo_covariates`**: A multi-layer `SpatRaster` with four environmental surfaces (`err27`, `ffp`, `cti`, `hli`) used in the `IBE_IBR_Framework` vignette.

### 13. Parallel processing improvements

The parallel evaluation of GA candidate solutions has been updated to correctly handle `SpatRaster` objects in parallel worker processes on Windows. Rasters are serialized to `PackedSpatRaster` (via `terra::wrap()`) before being sent to workers, then unwrapped inside each worker. The `parallel` argument in `GA.prep()` now also accepts a numeric value for an explicit core count.

### 14. Package renamed: ResistanceGA → ResistanceGA2

The package name changed from `ResistanceGA` to `ResistanceGA2`. Update all `library()` calls and internal references accordingly. The GitHub repository is `wpeterman/ResistanceGA2`.

---

## Quick Start

```r
library(ResistanceGA2)
library(terra)

# Package data
surfaces  <- resistance_surfaces          # 3-layer SpatRaster (unwrapped automatically)
pts_small <- vect(samples[, 2:3], type = "points")

# 1. Configure the GA
ga_inputs <- GA.prep(
  raster      = subset(surfaces, c("categorical", "continuous")),
  Results.dir = file.path(tempdir(), "ResistanceGA2-demo"),
  method      = "LL",
  monitor     = FALSE,
  quiet       = TRUE
)

# 2. Prepare distance engine (gdistance, no external software needed)
gdist_inputs <- gdist.prep(
  n.Pops  = nrow(samples),
  response = lower(as.matrix(dist(samples[, 2:3]))),
  samples = pts_small,
  method  = "costDistance"
)

# 3. Optimize single surfaces
ss_out <- SS_optim(
  gdist.inputs     = gdist_inputs,
  GA.inputs        = ga_inputs,
  dist_mod         = TRUE,
  null_mod         = TRUE,
  diagnostic_plots = FALSE
)

# 4. Examine results
summary(ss_out)
plot(ss_out)
ss_out$AICc
ss_out$ContinuousResults
```

For Julia/Circuitscape, replace `gdist.prep()` with `jl.prep()` and `gdist.inputs` with `jl.inputs`. See the `Julia_Guide` vignette.

---

## Distance Engines

| Engine | Function | Method | External dependency |
|--------|----------|--------|---------------------|
| gdistance | `gdist.prep()` | Least-cost path (`costDistance`) | None |
| gdistance | `gdist.prep()` | Random-walk commute (`commuteDistance`) | None |
| Julia/Circuitscape | `jl.prep()` | Circuit-theory resistance distance | Julia + Circuitscape.jl |

`commuteDistance` is the closest gdistance analog to Circuitscape's resistance distance. For iterative optimization over large landscapes or many sample populations, gdistance is typically faster. Use Julia/Circuitscape when you need Circuitscape-matched output or cumulative current maps.

---

## Vignettes

| Vignette | Description |
|----------|-------------|
| `GettingStarted` | Four-step workflow overview, package data, transformation inspection |
| `ResistanceGA` | Full tutorial: single-surface, multi-surface, all-combinations, and bootstrap analyses |
| `Surface_Transformations` | All eight transformation families, parameter guidance, and overfitting considerations |
| `IBE_IBR_Framework` | Joint IBD/IBE/IBR testing with the new `mlpe()` / `mlpe_data()` interface and Idaho demo data |
| `Julia_Guide` | Installing Julia, configuring JuliaConnectoR, running Circuitscape, troubleshooting |
| `Circuitscape_INI` | Direct Circuitscape runs with advanced INI file control |

```r
vignette("GettingStarted",        package = "ResistanceGA2")
vignette("ResistanceGA",          package = "ResistanceGA2")
vignette("Surface_Transformations", package = "ResistanceGA2")
vignette("IBE_IBR_Framework",     package = "ResistanceGA2")
vignette("Julia_Guide",           package = "ResistanceGA2")
```

---

## Notes for Existing ResistanceGA Users

- Replace `library(ResistanceGA)` with `library(ResistanceGA2)`.
- Replace `RasterLayer`/`RasterStack` with `SpatRaster` (`terra::rast()`).
- Replace `SpatialPoints` with `SpatVector` (`terra::vect()`).
- Replace `CS.prep()` / `Run_CS()` with `jl.prep()` / `Run_CS.jl()`.
- Remove any `k.smooth` arguments — this parameter no longer exists.
- Remove `Grid.Search()` calls — this function has been removed.
- `SS_optim.scale()` / `MS_optim.scale()` have been merged into `SS_optim()` / `MS_optim()`.
- Result objects now carry S3 classes; `class(ss_out)` returns a character vector rather than `"list"`.

---

## Recommended Practices

- **Use Monomolecular families by default.** For continuous surfaces, restrict `select.trans` to `"M"` unless there is a biological reason to expect a unimodal or quadratic resistance response. Ricker families increase overfitting risk.
- **Run replicates.** The GA is stochastic. Consistent parameter recovery across multiple replicates provides stronger support than a single run.
- **Report k.value.** Because parameter counting for pairwise MLPE models is not standardized, always report which `k.value` option was used. The default is `k.value = 2`.
- **Start small.** Test workflows with reduced `pop.size` (25) and `maxiter` (50–100) before committing to production settings.

---

## Citation

If you use this package, please cite:

Peterman, W. E. 2018. ResistanceGA: An R package for the optimization of resistance surfaces using genetic algorithms. *Methods in Ecology and Evolution* 9:1638–1647. <https://doi.org/10.1111/2041-210X.12984>

Additional relevant references:

Peterman, W. E., G. M. Connette, R. D. Semlitsch, and L. S. Eggert. 2014. Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. *Molecular Ecology* 23:2402–2413.

Ruiz-López, M. J., Barelli, C., Rovero, F., Hodges, K., Roos, C., Peterman, W. E., and Ting, N. 2016. A novel landscape genetic approach demonstrates the effects of human disturbance on the Udzungwa red colobus monkey (*Procolobus gordonorum*). *Heredity* 116:167–176.
