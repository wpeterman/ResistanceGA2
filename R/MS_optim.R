#' Simultaneous optimization of multiple resistance surfaces
#'
#' Optimize multiple resistance surfsaces simultaneously using genetic algorithms
#'
#' @param gdist.inputs Object from \code{\link{gdist.prep}}. Supply when optimizing with gdistance.
#' @param jl.inputs Object from \code{\link{jl.prep}}. Supply when optimizing with Circuitscape via Julia.
#' @param GA.inputs Object from \code{\link{GA.prep}}.
#' @param diagnostic_plots Logical. Generate and save diagnostic plots? Default = \code{TRUE}.
#' @return An object of class \code{resga_ms_optim}: a named list with GA
#'   summary, MLPE model, AICc table, distance matrices, percent contribution,
#'   and \code{k} table.
#' @usage MS_optim(gdist.inputs = NULL, jl.inputs = NULL, GA.inputs, diagnostic_plots = TRUE)

#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#' 
#' @examples
#' \dontrun{
#' pts <- terra::vect(sample_pops[[1]], type = "points")
#' gdist.inputs <- gdist.prep(
#'   n.Pops = nrow(sample_pops[[1]]),
#'   response = lower(Dc_list[[1]]),
#'   samples = pts
#' )
#'
#' GA.inputs <- GA.prep(
#'   raster = raster_orig,
#'   Results.dir = file.path(tempdir(), "ResistanceGA2-ms"),
#'   pop.size = 20,
#'   maxiter = 10,
#'   run = 5,
#'   quiet = TRUE,
#'   monitor = FALSE
#' )
#'
#' ms.out <- MS_optim(
#'   gdist.inputs = gdist.inputs,
#'   GA.inputs = GA.worker.inputs,
#'   diagnostic_plots = FALSE
#' )
#' }

MS_optim <- function(gdist.inputs = NULL,
                     jl.inputs    = NULL,
                     GA.inputs,
                     diagnostic_plots = TRUE) {
  k.value <- GA.inputs$k.value
  GA.worker.inputs <- .rga_prepare_parallel_inputs(GA.inputs)
  wd <- getwd()

  # gdistance ---------------------------------------------------------------
  
  if (!is.null(gdist.inputs)) {
    
    
    # MLPE with Covariates ----------------------------------------------------
    
    if(!is.null(gdist.inputs$covariates)) { 
      #  * Island GA -------------------------------------------------------------
      if(isTRUE(GA.inputs$gaisl)) {
        t1 <- proc.time()[3]
        multi.GA_nG <- gaisl(
          type = "real-valued",
          fitness = Resistance.Opt_multi.cov,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          mutation = GA.inputs$mutation,
          pcrossover = GA.inputs$pcrossover,
          crossover = GA.inputs$crossover,
          pmutation = GA.inputs$pmutation,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.worker.inputs,
          gdist.inputs = gdist.inputs,
          lower = GA.inputs$ga.min,
          upper = GA.inputs$ga.max,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          numIslands = GA.inputs$numIslands,
          migrationRate = GA.inputs$migrationRate,
          migrationInterval = GA.inputs$migrationInterval,
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          parallel = GA.inputs$parallel,
          run = GA.inputs$run,
          # keepBest = GA.inputs$keepBest,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          suggestions = GA.inputs$SUGGESTS,
          quiet = GA.inputs$quiet
        )
        rt <- proc.time()[3] - t1
        
      } else {
        # * Standard GA -------------------------------------------------------------
        
        t1 <- proc.time()[3]
        multi.GA_nG <- ga(
          type = "real-valued",
          fitness = Resistance.Opt_multi.cov,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          mutation = GA.inputs$mutation,
          pcrossover = GA.inputs$pcrossover,
          crossover = GA.inputs$crossover,
          pmutation = GA.inputs$pmutation,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.worker.inputs,
          gdist.inputs = gdist.inputs,
          lower = GA.inputs$ga.min,
          upper = GA.inputs$ga.max,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          parallel = GA.inputs$parallel,
          run = GA.inputs$run,
          keepBest = GA.inputs$keepBest,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          suggestions = GA.inputs$SUGGESTS,
          quiet = GA.inputs$quiet
        )
        rt <- proc.time()[3] - t1
        
      }
      
    } # End covariates
    
    # MLPE no Covariates ------------------------------------------------------
    if(is.null(gdist.inputs$covariates)) {
      #  * Island GA -------------------------------------------------------------
      if(isTRUE(GA.inputs$gaisl)) {
        t1 <- proc.time()[3]
        multi.GA_nG <- gaisl(
          type = "real-valued",
          fitness = Resistance.Opt_multi,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          mutation = GA.inputs$mutation,
          pcrossover = GA.inputs$pcrossover,
          crossover = GA.inputs$crossover,
          pmutation = GA.inputs$pmutation,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.worker.inputs,
          gdist.inputs = gdist.inputs,
          lower = GA.inputs$ga.min,
          upper = GA.inputs$ga.max,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          numIslands = GA.inputs$numIslands,
          migrationRate = GA.inputs$migrationRate,
          migrationInterval = GA.inputs$migrationInterval,
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          parallel = GA.inputs$parallel,
          run = GA.inputs$run,
          # keepBest = GA.inputs$keepBest,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          suggestions = GA.inputs$SUGGESTS,
          quiet = GA.inputs$quiet
        )
        rt <- proc.time()[3] - t1
        
      } else {
        # * Standard GA -------------------------------------------------------------
        
        t1 <- proc.time()[3]
        multi.GA_nG <- ga(
          type = "real-valued",
          fitness = Resistance.Opt_multi,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          mutation = GA.inputs$mutation,
          pcrossover = GA.inputs$pcrossover,
          crossover = GA.inputs$crossover,
          pmutation = GA.inputs$pmutation,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.worker.inputs,
          gdist.inputs = gdist.inputs,
          lower = GA.inputs$ga.min,
          upper = GA.inputs$ga.max,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          parallel = GA.inputs$parallel,
          run = GA.inputs$run,
          keepBest = GA.inputs$keepBest,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          suggestions = GA.inputs$SUGGESTS,
          quiet = GA.inputs$quiet
        )
        rt <- proc.time()[3] - t1
        
      }
    } # End no covariates
    
    
    
    multi.GA_nG.o <- multi.GA_nG
    
    saveRDS(multi.GA_nG, 
            file = paste0(GA.inputs$Results.dir, paste(GA.inputs$parm.type$name, collapse = "."), "_full.rds"))
    
    if(dim(multi.GA_nG@solution)[1] > 1) {
      multi.GA_nG@solution <- t(as.matrix(multi.GA_nG@solution[1,]))
    }
    
    Opt.parm <- GA.opt <- multi.GA_nG@solution
    
    for (i in 1:GA.inputs$n.layers) {
      if (GA.inputs$surface.type[i] == "cat") {
        ga.p <-
          GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        parm <- ga.p / min(ga.p)
        if(max(parm) > GA.inputs$max.cat){
          parm <- SCALE.vector(parm, 1, GA.inputs$max.cat)
        }
        Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                       1])] <- parm
        
      } else {
        parm <-
          GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                       1])] <- parm
      }
    }
    multi.GA_nG@solution <- Opt.parm
    
    RAST <-
      Combine_Surfaces(
        PARM = multi.GA_nG@solution,
        gdist.inputs = gdist.inputs,
        GA.inputs = GA.worker.inputs,
        rescale = TRUE,
        p.contribution = TRUE
      )
    
    p.cont <- RAST$percent.contribution
    RAST <- RAST$combined.surface
    
    NAME <- paste(GA.inputs$parm.type$name, collapse = ".")
    names(RAST) <- NAME
    cd <- Run_gdistance(gdist.inputs, RAST)
    dat <- gdist.inputs$df
    dat$cd <- scale(c(cd))
    
    write.table(
      as.matrix(cd),
      file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,"_distMat.csv"),
      sep = ",",
      row.names = F,
      col.names = F
    )
    terra::writeRaster(RAST,
                       paste0(GA.inputs$Results.dir, NAME, ".asc"),
                       overwrite = TRUE)

    type <- if (nrow(terra::unique(RAST)) > 15) "continuous" else "categorical"

    if(isTRUE(diagnostic_plots)){
      Diagnostic.Plots(
        resistance.mat = cd,
        genetic.dist = gdist.inputs$response,
        plot.dir = GA.inputs$Plots.dir,
        type = type,
        name = NAME,
        ID = gdist.inputs$ID,
        ZZ = gdist.inputs$ZZ
      )
    }
    
    if(!is.null(gdist.inputs$covariates)) { 
      MLPE.results <- NULL
    } else {
      # Get parameter estimates
      MLPE.results <- MLPE.lmm_coef(
        resistance = GA.inputs$Results.dir,
        genetic.dist = gdist.inputs$response,
        out.dir = GA.inputs$Results.dir,
        method = "gd",
        ID = gdist.inputs$ID,
        ZZ = gdist.inputs$ZZ
      )
    }
    
    fit.mod <- mlpe_rga(formula = gdist.inputs$formula,
                        data = dat,
                        ZZ = gdist.inputs$ZZ,
                        REML = FALSE)
    fit.mod_REML <- mlpe_rga(formula = gdist.inputs$formula,
                             data = dat,
                             ZZ = gdist.inputs$ZZ,
                             REML = TRUE)
    
    fit.stats <- suppressWarnings(r.squaredGLMM(
      fit.mod
    ))
    
    LL <- logLik(
      fit.mod
    )[[1]]
    
    MLPE.model <- fit.mod
    
    # MLPE.model <- MLPE.lmm(
    #   resistance = cd,
    #   pairwise.genetic = gdist.inputs$response,
    #   REML = F,
    #   ID = gdist.inputs$ID,
    #   ZZ = gdist.inputs$ZZ
    # )
    
    k <- .rga_multisurface_k(fit.mod, GA.inputs)
    
    n <- gdist.inputs$n.Pops
    aic <- (-2 * LL) + (2 * k)
    AICc <- 
      # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
      (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
    
    Result.txt(
      GA.results = multi.GA_nG,
      GA.inputs = GA.worker.inputs,
      method = gdist.inputs$method,
      Run.Time = rt,
      fit.stats = fit.stats,
      optim = GA.inputs$method,
      k = k,
      aic = aic,
      AICc = AICc,
      LL = LL[[1]],
      fit.mod_REML = fit.mod_REML
    )
    
    write.table(p.cont, file = paste0(GA.inputs$Results.dir, "Percent_Contribution.csv"), sep = ",",
                row.names = F,
                col.names = T)
    
    # save(multi.GA_nG, 
    #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
    
    saveRDS(multi.GA_nG, 
            file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
    
    # unlink(GA.inputs$Write.dir, recursive = T, force = T)
    
    k.df <- data.frame(
      surface = NAME,
      k = k,
      fixed.effects = .rga_fixed_effect_label(fit.mod)
    )
    
    cd.list <- list(as.matrix(cd))
    names(cd.list) <- NAME
    
    AICc.tab <- data.frame(surface = NAME,
                           fixed.effects = .rga_fixed_effect_label(fit.mod),
                           obj = multi.GA_nG@fitnessValue,
                           k = k,
                           AIC = aic,
                           AICc = AICc,
                           R2m = fit.stats[[1]],
                           R2c = fit.stats[[2]],
                           LL = LL)
    
    colnames(AICc.tab) <-
      c(
        "Surface",
        "Fixed.Effects",
        paste0("obj.func_", GA.inputs$method),
        "k",
        "AIC",
        "AICc",
        "R2m",
        "R2c",
        "LL"
      )
    
    out <- list(GA.summary = multi.GA_nG,
                MLPE.model = MLPE.model,
                MLPE.model_REML = fit.mod_REML,
                AICc.tab = AICc.tab,
                cd = cd.list,
                percent.contribution = p.cont,
                k = k.df)
    
    out <- resga_add_class(out, "resga_ms_optim")
    return(out)
  }
  
  
  # >>> Julia ---------------------------------------------------------------
  if (!is.null(jl.inputs)) {
    # setwd(jl.inputs$JULIA_HOME)
    
    # MLPE with Covariates ----------------------------------------------------
    
    if(!is.null(jl.inputs$covariates)) { 
      #  * Island GA -------------------------------------------------------------
      if(isTRUE(GA.inputs$gaisl)) {
        # stop("Optimization with covariates is not currently supported with gaisl!")
        
        t1 <- proc.time()[3]
        multi.GA_nG <- gaisl(
          type = "real-valued",
          fitness = Resistance.Opt_multi.cov,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          mutation = GA.inputs$mutation,
          pcrossover = GA.inputs$pcrossover,
          crossover = GA.inputs$crossover,
          pmutation = GA.inputs$pmutation,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.worker.inputs,
          jl.inputs = jl.inputs,
          lower = GA.inputs$ga.min,
          upper = GA.inputs$ga.max,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          numIslands = GA.inputs$numIslands,
          migrationRate = GA.inputs$migrationRate,
          migrationInterval = GA.inputs$migrationInterval,
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          parallel = GA.inputs$parallel,
          run = GA.inputs$run,
          # keepBest = GA.inputs$keepBest,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          suggestions = GA.inputs$SUGGESTS,
          quiet = GA.inputs$quiet
        )
        rt <- proc.time()[3] - t1
        
      } else {
        # * Standard GA -------------------------------------------------------------
        
        t1 <- proc.time()[3]
        multi.GA_nG <- ga(
          type = "real-valued",
          fitness = Resistance.Opt_multi.cov,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          mutation = GA.inputs$mutation,
          pcrossover = GA.inputs$pcrossover,
          crossover = GA.inputs$crossover,
          pmutation = GA.inputs$pmutation,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.worker.inputs,
          jl.inputs = jl.inputs,
          lower = GA.inputs$ga.min,
          upper = GA.inputs$ga.max,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          parallel = GA.inputs$parallel,
          run = GA.inputs$run,
          keepBest = GA.inputs$keepBest,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          suggestions = GA.inputs$SUGGESTS,
          quiet = GA.inputs$quiet
        )
        rt <- proc.time()[3] - t1
        
      }
    } # End covariates
    
    # MLPE no Covariates ------------------------------------------------------
    
    if (is.null(jl.inputs$covariates)) {
      #  * Island GA -------------------------------------------------------------
      if(isTRUE(GA.inputs$gaisl)) {
        t1 <- proc.time()[3]
        multi.GA_nG <- gaisl(
          type = "real-valued",
          fitness = Resistance.Opt_multi,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          mutation = GA.inputs$mutation,
          pcrossover = GA.inputs$pcrossover,
          crossover = GA.inputs$crossover,
          pmutation = GA.inputs$pmutation,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.worker.inputs,
          jl.inputs = jl.inputs,
          lower = GA.inputs$ga.min,
          upper = GA.inputs$ga.max,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          numIslands = GA.inputs$numIslands,
          migrationRate = GA.inputs$migrationRate,
          migrationInterval = GA.inputs$migrationInterval,
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          parallel = GA.inputs$parallel,
          run = GA.inputs$run,
          # keepBest = GA.inputs$keepBest,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          suggestions = GA.inputs$SUGGESTS,
          quiet = GA.inputs$quiet
        )
        rt <- proc.time()[3] - t1
        
      } else {
        # * Standard GA -------------------------------------------------------------
        
        t1 <- proc.time()[3]
        multi.GA_nG <- ga(
          type = "real-valued",
          fitness = Resistance.Opt_multi,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          mutation = GA.inputs$mutation,
          pcrossover = GA.inputs$pcrossover,
          crossover = GA.inputs$crossover,
          pmutation = GA.inputs$pmutation,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.worker.inputs,
          jl.inputs = jl.inputs,
          lower = GA.inputs$ga.min,
          upper = GA.inputs$ga.max,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          parallel = GA.inputs$parallel,
          run = GA.inputs$run,
          keepBest = GA.inputs$keepBest,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          suggestions = GA.inputs$SUGGESTS,
          quiet = GA.inputs$quiet
        )
        rt <- proc.time()[3] - t1
        
      }
    } # End no covariates
    
    
    multi.GA_nG.o <- multi.GA_nG
    
    saveRDS(multi.GA_nG, 
            file = paste0(GA.inputs$Results.dir, paste(GA.inputs$parm.type$name, collapse = "."), "_full.rds"))
    
    if(dim(multi.GA_nG@solution)[1] > 1) {
      multi.GA_nG@solution <- t(as.matrix(multi.GA_nG@solution[1,]))
    }
    
    Opt.parm <- GA.opt <- multi.GA_nG@solution
    for (i in 1:GA.inputs$n.layers) {
      if (GA.inputs$surface.type[i] == "cat") {
        ga.p <-
          GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        parm <- ga.p / min(ga.p)
        if(max(parm) > GA.inputs$max.cat){
          parm <- SCALE.vector(parm, 1, GA.inputs$max.cat)
        }
        Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                       1])] <- parm
        
      } else {
        parm <-
          GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                       1])] <- parm
      }
    }
    multi.GA_nG@solution <- Opt.parm
    
    RAST <-
      Combine_Surfaces(
        PARM = multi.GA_nG@solution,
        jl.inputs = jl.inputs,
        GA.inputs = GA.worker.inputs,
        rescale = TRUE,
        p.contribution = TRUE
      )
    
    p.cont <- RAST$percent.contribution
    RAST <- RAST$combined.surface
    
    NAME <- paste(GA.inputs$parm.type$name, collapse = ".")
    names(RAST) <- NAME
    
    cd <- suppressWarnings(Run_CS.jl(jl.inputs, RAST, full.mat = TRUE))
    cd.l <- lower(cd)
    cd.l <- cd.l[cd.l != -1]
    cd.l <- scale(cd.l)
    dat <- jl.inputs$df
    dat$cd <- cd.l
    # dat$cd <- scale(lower(cd)[which(lower(cd) != -1)])
    
    write.table(
      cd,
      file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
      sep = ",",
      row.names = F,
      col.names = F
    )
    terra::writeRaster(RAST,
                       paste0(GA.inputs$Results.dir, NAME, ".asc"),
                       overwrite = TRUE)

    type <- if (nrow(terra::unique(RAST)) > 15) "continuous" else "categorical"
    
    if(isTRUE(diagnostic_plots)){
      Diagnostic.Plots(
        resistance.mat = dat$cd,
        genetic.dist = jl.inputs$response,
        plot.dir = GA.inputs$Plots.dir,
        type = type,
        name = NAME,
        ID = jl.inputs$ID,
        ZZ = jl.inputs$ZZ
      )
    }
    
    # Get parameter estimates
    MLPE.results <- MLPE.lmm_coef(
      formula = jl.inputs$formula,
      inputs = dat,
      resistance = GA.inputs$Results.dir,
      genetic.dist = jl.inputs$response,
      out.dir = GA.inputs$Results.dir,
      method = "jl",
      ID = jl.inputs$ID,
      ZZ = jl.inputs$ZZ
    )
    
    # fit.stats <- r.squaredGLMM(
    #   MLPE.lmm(
    #     resistance = lower(cd),
    #     pairwise.genetic = jl.inputs$response,
    #     REML = F,
    #     ID = jl.inputs$ID,
    #     ZZ = jl.inputs$ZZ
    #   )
    # )
    # 
    # aic <- AIC(
    #   MLPE.lmm(
    #     resistance = lower(cd),
    #     pairwise.genetic = jl.inputs$response,
    #     REML = F,
    #     ID = jl.inputs$ID,
    #     ZZ = jl.inputs$ZZ
    #   )
    # )
    # 
    # LL <- logLik(
    #   MLPE.lmm(
    #     resistance = lower(cd),
    #     pairwise.genetic = jl.inputs$response,
    #     REML = F,
    #     ID = jl.inputs$ID,
    #     ZZ = jl.inputs$ZZ
    #   )
    # )
    # 
    # MLPE.model <- MLPE.lmm(
    #   resistance = lower(cd),
    #   pairwise.genetic = jl.inputs$response,
    #   REML = F,
    #   ID = jl.inputs$ID,
    #   ZZ = jl.inputs$ZZ
    # )
    
    fit.mod <- mlpe_rga(formula = jl.inputs$formula,
                        data = dat,
                        ZZ = jl.inputs$ZZ,
                        REML = FALSE)
    fit.mod_REML <- mlpe_rga(formula = jl.inputs$formula,
                             data = dat,
                             ZZ = jl.inputs$ZZ,
                             REML = TRUE)
    
    fit.stats <- suppressWarnings(r.squaredGLMM(
      fit.mod
    ))
    
    LL <- logLik(
      fit.mod
    )[[1]]
    
    MLPE.model <- fit.mod
    
    k <- .rga_multisurface_k(fit.mod, GA.inputs)
    
    n <- jl.inputs$n.Pops
    aic <- (-2 * LL) + (2 * k)
    AICc <- 
      # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
      (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
    
    Result.txt(
      GA.results = multi.GA_nG,
      GA.inputs = GA.worker.inputs,
      method = "CIRCUITSCAPE.jl",
      Run.Time = rt,
      fit.stats = fit.stats,
      optim = GA.inputs$method,
      k = k,
      aic = aic,
      AICc = AICc,
      LL = LL[[1]],
      fit.mod_REML = fit.mod_REML
    )
    
    write.table(p.cont, file = paste0(GA.inputs$Results.dir, "Percent_Contribution.csv"), sep = ",",
                row.names = F,
                col.names = T)
    
    # save(multi.GA_nG, 
    #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
    
    saveRDS(multi.GA_nG, 
            file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
    
    # unlink(GA.inputs$Write.dir, recursive = T, force = T)
    
    k.df <- data.frame(
      surface = NAME,
      k = k,
      fixed.effects = .rga_fixed_effect_label(fit.mod)
    )
    
    cd.list <- list(as.matrix(cd))
    names(cd.list) <- NAME
    
    AICc.tab <- data.frame(surface = NAME,
                           fixed.effects = .rga_fixed_effect_label(fit.mod),
                           obj = multi.GA_nG@fitnessValue,
                           k = k,
                           AIC = aic,
                           AICc = AICc,
                           R2m = fit.stats[[1]],
                           R2c = fit.stats[[2]],
                           LL = LL)
    
    colnames(AICc.tab) <-
      c(
        "Surface",
        "Fixed.Effects",
        paste0("obj.func_", GA.inputs$method),
        "k",
        "AIC",
        "AICc",
        "R2m",
        "R2c",
        "LL"
      )
    
    setwd(wd)
    
    out <- list(GA.summary = multi.GA_nG,
                MLPE.model = MLPE.model,
                MLPE.model_REML = fit.mod_REML,
                AICc.tab = AICc.tab,
                cd = cd.list,
                percent.contribution = p.cont,
                k = k.df)
    
    out <- resga_add_class(out, "resga_ms_optim")
    return(out)
  } # End Julia
} # End function
