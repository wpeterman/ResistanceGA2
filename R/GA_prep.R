#' Create R object with genetic algorithm optimization settings
#'
#' This function prepares and compiles objects and commands for optimization with the GA package
#'
#' @param raster A \code{SpatRaster} (from \pkg{terra}) with one or more layers.
#' @param Results.dir Directory to export optimization results. If the directory
#'   does not exist, it will be created automatically. If it exists and is not
#'   empty, a message will note that existing results may be overwritten. It is
#'   critical that there are NO SPACES in the directory path. If using the
#'   \code{\link[ResistanceGA2]{all_comb}} function, specify \code{Results.dir}
#'   as \code{"all_comb"}.
#' @param min.cat The minimum value to be assessed during optimization of categorical resistance surfaces (Default = 1 / max.cat)
#' @param max.cat The maximum value to be assessed during optimization of categorical resistance surfaces (Default = 1000)
#' @param max.cont The maximum value to be assessed during optimization of continuous resistance surfaces (Default = 1000)
#' @param shape.min The minimum value for the shape parameter used for transforming resistance surfaces. If unspecified, used 0.5
#' @param shape.max The maximum value for the shape parameter used for transforming resistance surfaces. If unspecified, used 14.5
#' @param cont.shape A vector of hypothesized relationships that each continuous resistance surface will have in relation to the genetic distance response (Default = NULL; see details)
#' @param select.trans Option to specify which transformations are applied to
#'   continuous surfaces. Must be provided as a list. "A" = All, "M" =
#'   Monomolecular only, "R" = Ricker only. Default = "M"; see Details.
#' @param cat.levels Number of unique levels to permit in categorical surface (Default = 15). See Details
#' @param method Objective function to be optimized. Select "AIC", "R2", or "LL" to optimize resistance surfaces based on AIC, variance explained (R2), or log-likelihood. (Default = "LL")
#' @param k.value Specification of how k, the number of parameters in the mixed effects model, is determined. Specify 1, 2, 3, or 4 (Default = 2; see details).
#'
#' 1 --> k = 2;
#'
#' 2 --> k = number of parameters optimized plus the intercept;
#'
#' 3 --> k =  the number of parameters optimized plus the intercept and the number of layers optimized;
#' 
#' 4 --> k = the number of layers optimized plus the intercept
#' @param pop.mult Value will be multiplied with number of parameters in surface to determine 'popSize' in GA. By default this is set to 15.
#' @param percent.elite An integer percent used to determine the number of best fitness individuals to survive at each generation ('elitism' in GA). By default the top 5\% individuals will survive at each iteration.
#' @param type Default is "real-valued"
#' @param population Default is "gareal_Population" from GA
#' @param selection Default is "gareal_lsSelection" from GA
#' @param mutation Default is "gareal_raMutation" from GA
#' @param pcrossover Probability of crossover. Default = 0.85
#' @param pmutation Probability of mutation. Default = 0.125
#' @param crossover Default = "gareal_blxCrossover". This crossover method greatly improved optimization during preliminary testing
#' @param maxiter Maximum number of iterations to run before the GA search is halted. If using standard \code{ga} optimizer, the default = 1000. If using \code{gaisl = TRUE}, then this is set to 15x the \code{migrationInterval}
#' @param pop.size Number of individuals to create each generation. If \code{gaisl = TRUE}, then this number is automatically calculated as \code{numIslands} * \code{island.pop} 
#' @param parallel A logical argument specifying if parallel computing should be used (TRUE) or not (FALSE, default) for evaluating the fitness function. You can also specify the number of cores to use. 
#' @param gaisl Should the genetic algorithm use the islands parallel optimization? (Default = FALSE)
#' @param island.pop The number of individuals to populate each island. (Default = 20)
#' @param numIslands If \code{gaisl = TRUE}, an integer value which specifies the number of islands to use in the genetic evolution (by default will be set to 4)
#' @param migrationRate If \code{gaisl = TRUE}, a value in the range (0, 1) which gives the proportion of individuals that undergo migration between islands in every exchange (by default equal to 0.10).
#' @param migrationInterval If \code{gaisl = TRUE}, an integer value specifying the number of iterations at which exchange of individuals takes place. This interval between migrations is called an epoch, and it is set at 10 by default.
#' @param run Number of consecutive generations or epochs without any improvement in objective function before the GA is stopped. If using standard \code{ga}, the default = 25. If using \code{gaisl = TRUE}, then the default \code{run} value will be calculated as \code{migrationInterval} * 5
#' @param keepBest A logical argument specifying if best solutions at each iteration should be saved (Default = TRUE)
#' @param optim A logical defaulting to \code{FALSE} determining whether or not a local search using general-purpose optimisation algorithms should be used. See argument \code{optimArgs} for further details and finer control. Setting to TRUE has the potential to improve optimization accuracy, but will increase optimization time.
#' @param optim.method The method to be used among those available in \code{\link[stats]{optim}} function. By default, the BFGS algorithm with box constraints is used, where the bounds are those provided in the \code{ga()} function call. Further methods are available as described in the Details section in help(optim).
#' @param poptim A value in the range [0,1] specifying the probability of performing a local search at each iteration of GA (default 0.0). Only change if your optimization is relatively fast.
#' @param pressel A value in the range [0,1] specifying the pressure selection (default 1.00). The local search is started from a random solution selected with probability proportional to fitness. High values of pressel tend to select the solutions with the largest fitness, whereas low values of pressel assign quasi-uniform probabilities to any solution.
#' @param control A list of control parameters. See 'Details' section in \code{\link[stats]{optim}}
#' @param hessian	Logical. Should a numerically differentiated Hessian matrix be returned? This will allow for the calculation of standard errors on parameter estimates (not yet implemented). Default = FALSE
#' @param opt.digits The number of significant digits that the objective function will be assessed at. By default, no rounding occurs.
#' @param seed Integer random number seed to replicate \code{ga} optimization
#' @param monitor Default = TRUE, which prints the average and best fitness values at each iteration. 
#' @param quiet Logical. If TRUE, the objective function and step run time will not be printed to the screen after each step. Only \code{ga} summary information will be printed following each iteration. (Default = FALSE)
#' @return An R object that is a required input into optimization functions
#'
#' @details Provide the resistance surfaces directly as a \code{SpatRaster}. If
#' multiple surfaces are to be optimized, include them as layers in that
#' object. File-path and directory-based raster inputs are no longer accepted.
#'
#' The Default for \code{k.value} is 2, which sets k equal to the number of parameters optimized, plus 1 for the intercept term. Prior to version 3.0-0, \code{k.value} could not be specified by the user and followed setting 2, such that k was equal to the number of parameters optimized plus the intercept term.
#'
#' \code{cont.shape} can take values of "Increase", "Decrease", or "Peaked". If you believe a resistance surface is related to your response in a particular way, specifying this here may decrease the time to optimization. \code{cont.shape} is used to generate an initial set of parameter values to test during optimization. If specified, a greater proportion of the starting values will include your believed relationship. If unspecified (the Default), a completely random set of starting values will be generated.
#'
#' If it is desired that only certain transformations be assessed for
#' continuous surfaces, then this can be specified using
#' \code{select.trans}. By default, only monomolecular transformations
#' will be assessed for continuous surfaces unless otherwise specified.
#' This default is also the recommended starting point for most analyses:
#' monomolecular families are usually easier to interpret and generally
#' reduce the risk of overfitting relative to Ricker families. Ricker
#' transformations should usually be explored only when there is a strong
#' biological or ecological justification for a unimodal or quadratic-type
#' relationship. Specific transformations can be specified by providing a
#' vector of values (e.g., \code{c(1,3,5)}), with values corresponding to the
#' equation numbers as detailed in \code{\link[ResistanceGA2]{Resistance.tran}}.
#' If multiple rasters are to be optimized from the same \code{SpatRaster}
#' object, then a list of transformations must be provided in the order that
#' the raster surfaces will be assessed. For example:\cr
#' \code{select.trans = list("M", "A", "R", c(5,6))}\cr
#' will result in surface one only being optimized with Monomolecular transformations, surface two with all transformations, surface three with only Ricker transformations, and surface four with Reverse Ricker and Reverse Monomolecular only. If a categorical surface is among the rasters to be optimized, it is necessary to specify \code{NA} to accommodate this.
#' 
#' \code{cat.levels} defaults to 15. This means that when a raster surface has <= 15 unique levels, it will be treated as a categorical surface in the analysis. This value can be increased, but optimization of surfaces with many levels may take more time. Additionally, depending upon the prevalence and configuration of categorical features and spatial sample locations, some levels are likely to be poorly estimated. This may be evident if estimated resistance values vary substantially between runs of ResistanceGA.
#' 
#' Setting \code{gaisl = TRUE} has the potential greatly reduce the optimization run time, potentially with greater accuracy. This is a distributed multiple-population GA, where the population is partitioned into several subpopulations and assigned to separated islands. Independent GAs are executed in each island, and only occasionally sparse exchanges of individuals are performed among the islands. 
#'
#' It is recommended to first run GA optimization with the default settings

#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#' @usage GA.prep(raster,
#'                Results.dir = NULL,
#'                min.cat = NULL,
#'                max.cat = 1000,
#'                max.cont = 1000,
#'                shape.min = NULL,
#'                shape.max = NULL,
#'                cont.shape = NULL,
#'                select.trans = NULL,
#'                cat.levels = 15,
#'                method = "LL",
#'                k.value = 2,
#'                pop.mult = 15,
#'                percent.elite = 5,
#'                type = "real-valued",
#'                pcrossover = 0.85,
#'                pmutation = 0.125,
#'                maxiter = 1000,
#'                run = NULL,
#'                keepBest = TRUE,
#'                population = gaControl(type)$population,
#'                selection = gaControl(type)$selection,
#'                crossover = "gareal_blxCrossover",
#'                mutation = gaControl(type)$mutation,
#'                pop.size = NULL,
#'                parallel = FALSE,
#'                gaisl = FALSE,
#'                island.pop = 20,
#'                numIslands = NULL,
#'                migrationRate = NULL,
#'                migrationInterval = NULL,
#'                optim = FALSE,
#'                optim.method = "L-BFGS-B", 
#'                poptim = 0.0,
#'                pressel = 1.00,
#'                control = list(fnscale = -1, maxit = 100),
#'                hessian = FALSE,
#'                opt.digits = NULL,
#'                seed = NULL,
#'                monitor = TRUE,
#'                quiet = FALSE)
#' 
#' @examples
#' ga.inputs <- GA.prep(
#'   raster = raster_orig,
#'   Results.dir = file.path(tempdir(), "ResistanceGA2-example"),
#'   quiet = TRUE,
#'   monitor = FALSE
#' )
#'
#' ga.inputs$layer.names

GA.prep <- function(raster,
                    Results.dir = NULL,
                    min.cat = NULL,
                    max.cat = 1000,
                    max.cont = 1000,
                    shape.min = NULL,
                    shape.max = NULL,
                    cont.shape = NULL,
                    select.trans = NULL,
                    cat.levels = 15,
                    method = "LL",
                    k.value = 2,
                    pop.mult = 15,
                    percent.elite = 5,
                    type = "real-valued",
                    pcrossover = 0.85,
                    pmutation = 0.125,
                    maxiter = 1000,
                    run = NULL,
                    keepBest = TRUE,
                    population = gaControl(type)$population,
                    selection = gaControl(type)$selection,
                    crossover = "gareal_blxCrossover",
                    mutation = gaControl(type)$mutation,
                    pop.size = NULL,
                    parallel = FALSE,
                    gaisl = FALSE,
                    island.pop = 20,
                    numIslands = NULL,
                    migrationRate = NULL,
                    migrationInterval = NULL,
                    optim = FALSE,
                    optim.method = "L-BFGS-B", 
                    poptim = 0.0,
                    pressel = 1.00,
                    control = list(fnscale = -1, maxit = 100),
                    hessian = FALSE,
                    opt.digits = NULL,
                    seed = NULL,
                    monitor = TRUE,
                    quiet = FALSE) {
  
  
  
  # Hybrid GA ---------------------------------------------------------------
  if(isTRUE(optim)) {
    optimArgs <- list(method = optim.method,
                      poptim = poptim,
                      pressel = pressel,
                      control = control,
                      hessian = hessian)
  } else {
    optimArgs <- NULL
  }
  
  
  # Island model ------------------------------------------------------------
  
  if(isTRUE(gaisl)) {
    if(is.numeric(parallel)) {
      parallel <- parallel
    } else {
      parallel <- TRUE
    }
    
    if(is.null(numIslands)) {
      numIslands <- 4
    }
    
    if(is.null(migrationRate)) {
      migrationRate <- 0.10
    }
    
    if(is.null(migrationInterval)) {
      migrationInterval <- 10
    }
    
    if(is.null(maxiter)) {
      maxiter <- migrationInterval * 15
    }
    
    if(is.null(pop.size)) {
      pop.size <- island.pop * numIslands
    }
    
    if(is.null(run)) {
      run <- migrationInterval * 5
    }
  } # End island set up
  
  if(is.null(run)) {
    run <- 25
  }
  
  if(is.null(min.cat)){
    min.cat <- 1 / max.cat
  }
  
  if(is.null(shape.min)){
    shape.min <- 0.5
  }
  
  if(is.null(shape.max)){
    shape.max <- 14.5
  }
  
  if(is.null(Results.dir)) {
    warning(paste0(
      "'Results.dir' was not specified. Results will be exported to ",
      getwd()
    ))
    Results.dir <- paste0(getwd(), "/")
  }

  if ((Results.dir != 'all.comb') & (Results.dir != 'all_comb')) {
    if (!dir.exists(Results.dir)) {
      dir.create(Results.dir, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(Results.dir)) {
        stop("Failed to create 'Results.dir': ", Results.dir)
      }
      message("Created 'Results.dir': ", Results.dir)
    } else {
      existing <- list.files(Results.dir, all.files = FALSE, no.. = TRUE)
      if (length(existing) > 0) {
        message("NOTE: 'Results.dir' is not empty. Existing results may be overwritten.")
      }
    }

    results_root <- paste0(
      sub("[/\\\\]+$", "", normalizePath(Results.dir, winslash = "/", mustWork = FALSE)),
      "/"
    )
  }
  
  if(!is.null(select.trans)){
    if(!is.list(select.trans)){
      stop("Select transformations must by provided as a list. See Details.")
    }
  }
  
  
  r <- .validate_spatraster(raster, arg = "raster")
  layer.names <- names(r)
  n.layers <- terra::nlyr(r)
  
  if(Results.dir != 'all.comb' & Results.dir != 'all_comb') {
    if ("Results" %in% dir(Results.dir) == FALSE)
      dir.create(file.path(Results.dir, "Results"))
    Results.DIR <- paste0(results_root, "Results/")
    
    # if ("tmp" %in% dir(Results.dir) == FALSE)
    #   dir.create(file.path(Results.dir, "tmp"))
    # Write.dir <- paste0(Results.dir, "tmp/")
    Write.dir <- paste0(normalizePath(tempdir()),'/')
    
    
    if ("Plots" %in% dir(Results.DIR) == FALSE)
      dir.create(file.path(Results.DIR, "Plots"))
    Plots.dir <- paste0(Results.DIR, "Plots/")
  }
  
  if(Results.dir == 'all.comb' | Results.dir == 'all_comb') {
    Results.DIR <- NULL
    Write.dir <- NULL
    Plots.dir <- NULL
  }
  
  
  # Determine total number of parameters and types of surfaces included
  parm.type <- data.frame()
  min.list <- list()
  max.list <- list()
  SUGGESTS <- list()
  eqs <- list()
  for (i in 1:n.layers) {
    n.levels <- nrow(terra::unique(r[[i]]))
    
    if (n.levels <= cat.levels) {
      Level.val <- terra::unique(r[[i]])[[1]]
      parm.type[i, 1] <- "cat"
      parm.type[i, 2] <- n.levels
      parm.type[i, 3] <- layer.names[i]
      min.list[[i]] <- c(1, rep(min.cat, (n.levels - 1)))
      max.list[[i]] <- c(1, rep(max.cat, (n.levels - 1)))
      
      eqs[[i]] <- NA
      
    } else {
      parm.type[i, 1] <- "cont"
      parm.type[i, 2] <- 3
      min.list[[i]] <- c(1, shape.min, 0.001)  # eq, shape, max
      max.list[[i]] <- c(9.99, shape.max, max.cont)

      parm.type[i, 3] <- layer.names[i]
      
      if (is.null(select.trans)) {
        eqs[[i]] <- eq.set("M")
      } else {
        eqs[[i]] <- eq.set(select.trans[[i]])
      }
    }
  }
  
  
  colnames(parm.type) <- c("type", "n.parm", "name")
  parm.index <- c(0, cumsum(parm.type$n.parm))
  ga.min <- unlist(min.list)
  ga.max <- unlist(max.list)
  surface.type <- parm.type$type
  
  if (is.null(pop.size)) {
    if (length(ga.min) < 10) {
      pop.size <- min(c(15 * length(ga.min), 100))
    } else {
      pop.size <- 10 * length(ga.min)
    }
  }
  
  for (i in 1:length(surface.type)) {
    if (surface.type[i] == "cat") {
      SUGGESTS[[i]] <- sv.cat(levels = parm.type[i, 2],
                              pop.size = pop.size,
                              min.cat,
                              max.cat)
    } else if (exists("cont.shape") && length(cont.shape) > 0) {
      SUGGESTS[[i]] <- sv.cont.nG(cont.shape[1], pop.size = pop.size, max.cont, eqs = eqs[[i]])
      cont.shape <- cont.shape[-1]
    } else {
      SUGGESTS[[i]] <- sv.cont.nG("NA", pop.size = pop.size, max.cont, eqs = eqs[[i]])
    }
  }
  SUGGESTS <-
    matrix(unlist(SUGGESTS),
           nrow = nrow(SUGGESTS[[1]]),
           byrow = F)
  
  if (method != "AIC") {
    Min.Max <- 'min'
  } else {
    Min.Max <- 'max'
  }
  
  if(Results.dir != "all.comb" & Results.dir != "all_comb") {
    list(
      parm.index = parm.index,
      ga.min = ga.min,
      ga.max = ga.max,
      select.trans = eqs,
      surface.type = surface.type,
      parm.type = parm.type,
      Resistance.stack = r,
      n.layers = n.layers,
      layer.names = layer.names,
      pop.size = pop.size,
      min.list = min.list,
      max.list = max.list,
      SUGGESTS = SUGGESTS,
      raster = r,
      Results.dir = Results.DIR,
      Write.dir = Write.dir,
      Plots.dir = Plots.dir,
      type = type,
      pcrossover = pcrossover,
      pmutation = pmutation,
      crossover = crossover,
      maxiter = maxiter,
      run = run,
      keepBest = keepBest,
      population = population,
      selection = selection,
      mutation = mutation,
      parallel = parallel,
      gaisl = gaisl,
      numIslands = numIslands,
      migrationRate = migrationRate,
      migrationInterval = migrationInterval,
      optim = optim,
      optimArgs = optimArgs,
      pop.mult = pop.mult,
      percent.elite = percent.elite,
      Min.Max = Min.Max,
      max.cat = max.cat,
      method = method,
      k.value = k.value,
      opt.digits = opt.digits,
      seed = seed,
      monitor = monitor,
      quiet = quiet
    )
  } else {
    list(
      parm.index = parm.index,
      ga.min = ga.min,
      ga.max = ga.max,
      select.trans = eqs,
      surface.type = surface.type,
      parm.type = parm.type,
      Resistance.stack = r,
      n.layers = n.layers,
      layer.names = layer.names,
      pop.size = pop.size,
      min.list = min.list,
      max.list = max.list,
      SUGGESTS = SUGGESTS,
      raster = r,
      Results.dir = Results.DIR,
      Write.dir = Write.dir,
      Plots.dir = Plots.dir,
      type = type,
      pcrossover = pcrossover,
      pmutation = pmutation,
      crossover = crossover,
      maxiter = maxiter,
      run = run,
      keepBest = keepBest,
      population = population,
      selection = selection,
      mutation = mutation,
      parallel = parallel,
      gaisl = gaisl,
      numIslands = numIslands,
      migrationRate = migrationRate,
      migrationInterval = migrationInterval,
      optim = optim,
      optimArgs = optimArgs,
      pop.mult = pop.mult,
      percent.elite = percent.elite,
      Min.Max = Min.Max,
      max.cat = max.cat,
      method = method,
      k.value = k.value,
      opt.digits = opt.digits,
      seed = seed,
      monitor = monitor,
      quiet = quiet,
      inputs = list(
        raster = r,
        Results.dir = Results.dir,
        min.cat = min.cat,
        max.cat = max.cat,
        max.cont = max.cont,
        cont.shape = cont.shape,
        select.trans = select.trans,
        method = method,
        k.value = k.value,
        pop.mult = pop.mult,
        percent.elite = percent.elite,
        type = type,
        pcrossover = pcrossover,
        pmutation = pmutation,
        maxiter = maxiter,
        run = run,
        keepBest = keepBest,
        population = population,
        selection = selection,
        crossover = crossover,
        mutation = mutation,
        pop.size = pop.size,
        parallel = parallel,
        gaisl = gaisl,
        island.pop = island.pop,
        numIslands = numIslands,
        migrationRate = migrationRate,
        migrationInterval = migrationInterval,
        optim = optim,
        optim.method = optim.method,
        poptim = poptim,
        pressel = pressel,
        control = control,
        hessian = hessian,
        opt.digits = opt.digits,
        seed = seed,
        monitor = monitor,
        quiet = quiet
      )
    )
  }
  
}
