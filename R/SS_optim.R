#' Single surface optimization
#'
#' Optimize each surface in a \code{SpatRaster} stack using a genetic algorithm
#' executed with the \code{\link[GA]{ga}} function in \pkg{GA}.
#'
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA2]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param jl.inputs Object created from running \code{\link[ResistanceGA2]{jl.prep}} function. Defined if optimizing using CIRCUITSCAPE run in Julia
#' @param GA.inputs Object created from running \code{\link[ResistanceGA2]{GA.prep}} function
#' @param dist_mod Logical, if TRUE, a Distance model will be calculated and added to the output table (default = TRUE)
#' @param null_mod Logical, if TRUE, an intercept-only model will be calculated and added to the output table (default = TRUE)
#' @param diagnostic_plots Plotting and saving of diagnostic plots (Default = TRUE)
#' @return An object of class \code{resga_ss_optim}. This function optimizes
#' resistance surfaces in isolation. Following optimization of all surfaces,
#' several summary objects are created.\cr
#' \enumerate{
#' \item Diagnostic plots of model fit are output to the \code{Results/Plots}
#' directory created beneath \code{Results.dir}.
#' \item A .csv file with the Maximum Likelihood Population Effects mixed effects model coefficient estimates (MLPE_coeff_Table.csv)
#' \item Three summary .csv files are generated: CategoricalResults.csv, ContinuousResults.csv, & All_Results_AICc.csv. These tables contain AICc values and optimization summaries for each surface.
#' }
#' All results tables are also summarized in a named list ($ContinuousResults, $CategoricalResults, $AICc, $MLPE, $MLPE.list, $cd, $k)\cr
#' The \code{lmer} model objects stored $MLPE.list are fit using Restricted Maximum Likelihood \cr
#' $cd is a list of the optimized cost pairwise distance matrices and $k is a table of the surface names and number of parameters used to calculate AICc. These two objects can be passed to \code{\link[ResistanceGA2]{Resist.boot}} to conduct a bootstrap analysis.
#' @usage SS_optim(gdist.inputs,
#'  jl.inputs,
#'  GA.inputs,
#'  dist_mod, 
#'  null_mod,
#'  diagnostic_plots = TRUE)
#' @author Bill Peterman <Peterman.73@@osu.edu>
#' @export
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
#'   Results.dir = file.path(tempdir(), "ResistanceGA2-ss"),
#'   pop.size = 20,
#'   maxiter = 10,
#'   run = 5,
#'   quiet = TRUE,
#'   monitor = FALSE
#' )
#'
#' ss.out <- SS_optim(
#'   gdist.inputs = gdist.inputs,
#'   GA.inputs = GA.inputs,
#'   dist_mod = FALSE,
#'   null_mod = FALSE,
#'   diagnostic_plots = FALSE
#' )
#' }

SS_optim <- function(gdist.inputs = NULL,
                     jl.inputs = NULL,
                     GA.inputs,
                     dist_mod = TRUE,
                     null_mod = TRUE,
                     diagnostic_plots = TRUE) {
  
  t1 <- proc.time()[3]
  RESULTS.cat <- list() # List to store categorical results within
  RESULTS.cont <- list() # List to store continuous results within
  cnt1 <- 0
  cnt2 <- 0
  k.value <- GA.inputs$k.value
  GA.worker.inputs <- .rga_prepare_parallel_inputs(GA.inputs)
  MLPE.list <- list()
  cd.list <- list()
  k.list <- list()
  ga.list <- list()
  
  wd <- getwd()
  
  # Optimize each surface in turn
  for (i in 1:GA.inputs$n.layers) {
    r <- GA.inputs$Resistance.stack[[i]]
    names(r) <- GA.inputs$layer.names[i]
    
    
    # >>> gdistance <<< -------------------------------------------------
    
    if (!is.null(gdist.inputs)) {
      
      # MLPE with Covariates ----------------------------------------------------
      
      if(!is.null(gdist.inputs$covariates)) { 
        
        # * Island GA ---------------------------------------------------------------
        
        
        if(isTRUE(GA.inputs$gaisl)) {
          
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              SCALE.vector(single.GA@solution, 1, GA.inputs$max.cat)
            
            # single.GA@solution <-
            #   single.GA@solution / min(single.GA@solution)
            df <- data.frame(id = terra::unique(r)[[1]], t(single.GA@solution))
            r <- terra::subst(r, from = df$id, to = df[[2]])
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,  "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = cd,
                genetic.dist = gdist.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "categorical",
                name = NAME,
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
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 
                length(lme4::fixef(fit.mod)) - 1
              
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + 
                length(GA.inputs$layer.names) + 
                length(lme4::fixef(fit.mod)) - 1
              
            } else {
              k <- length(GA.inputs$layer.names[i]) + 
                length(lme4::fixef(fit.mod)) - 1
              
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <- fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # rm(single.GA, r)
            gc()
            
            # *-* Continuous -----------------------------------------------------------
            
          } else {
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method, "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = cd,
                genetic.dist = gdist.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "continuous",
                name = NAME,
                ID = gdist.inputs$ID,
                ZZ = gdist.inputs$ZZ
              )
              
              Plot.trans(
                PARM = single.GA@solution[-1],
                Resistance = GA.inputs$Resistance.stack[[i]],
                transformation = EQ,
                print.dir = GA.inputs$Plots.dir
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
            
            
            fit.stats <- r.squaredGLMM(
              fit.mod
            )
            
            LL <- logLik(
              fit.mod
            )[[1]]
            
            MLPE.list[[i + 1]] <- fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 
                length(lme4::fixef(fit.mod)) - 1
              
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + 
                length(GA.inputs$layer.names) + 
                length(lme4::fixef(fit.mod)) - 1
              
            } else {
              k <- length(GA.inputs$layer.names[i]) + 
                length(lme4::fixef(fit.mod)) - 1
              
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
            # rm(single.GA, r)
            gc()
          } # Close gaisl cat-cont if else
        } else { # * Standard GA -------------------------------------------------------------
          
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              parallel = GA.inputs$parallel,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              SCALE.vector(single.GA@solution, 1, GA.inputs$max.cat)
            
            # single.GA@solution <-
            #   single.GA@solution / min(single.GA@solution)
            # df <- data.frame(id = terra::unique(r)[[1]], t(single.GA@solution))
            # r <- subs(r, df)
            df <- data.frame(id = terra::unique(GA.inputs$Resistance.stack[[i]])[[1]], t(single.GA@solution))
            
            r <- terra::subst(GA.inputs$Resistance.stack[[i]], from = df$id, to = df[[2]]) ## Modified 2019-03-26
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,  "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA, 
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = cd,
                genetic.dist = gdist.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "categorical",
                name = NAME,
                ID = gdist.inputs$ID,
                ZZ = gdist.inputs$ZZ
              )
            }
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # if (k.value == 1) {
            #   k <- 2
            # } else if (k.value == 2) {
            #   k <- GA.inputs$parm.type$n.parm[i] + 1
            # } else if (k.value == 3) {
            #   k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            # } else {
            #   k <- length(GA.inputs$layer.names[i]) + 1
            # }
            
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
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + length(lme4::fixef(fit.mod)) - 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + length(lme4::fixef(fit.mod)) - 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + length(lme4::fixef(fit.mod)) - 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # rm(single.GA, r)
            gc()
          } 
          else { # *-* Continuous ----------
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method, "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA, 
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = cd,
                genetic.dist = gdist.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "continuous",
                name = NAME,
                ID = gdist.inputs$ID,
                ZZ = gdist.inputs$ZZ
              )
              
              Plot.trans(
                PARM = single.GA@solution[-1],
                Resistance = GA.inputs$Resistance.stack[[i]],
                transformation = EQ,
                print.dir = GA.inputs$Plots.dir
              )
            }
            
            
            # fit.stats <-
            #   r.squaredGLMM(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            # 
            # 
            # aic <-
            #   AIC(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            # 
            # LL <-
            #   logLik(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            
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
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
            # rm(single.GA, r)
            gc()
          } # Close cat-cont if else
          
        } # Close GA vs. gaisl if-else
        
        if (dist_mod == TRUE) {
          rd <- terra::classify(r, matrix(c(-Inf, Inf, 1), ncol = 3))
          names(rd) <- "dist"
          cd <- Run_gdistance(gdist.inputs, rd)
          
          dat <- gdist.inputs$df
          dat$cd <- scale(c(cd))
          
          write.table(
            as.matrix(cd),
            file = paste0(GA.inputs$Results.dir, 'Distance', "_", gdist.inputs$method, "_distMat.csv"),
            sep = ",",
            row.names = F,
            col.names = F
          )
          
          fit.mod <- mlpe_rga(formula = gdist.inputs$formula,
                              data = dat,
                              ZZ = gdist.inputs$ZZ,
                              REML = FALSE)
          
          fit.mod_REML <- mlpe_rga(formula = gdist.inputs$formula,
                                   data = dat,
                                   ZZ = gdist.inputs$ZZ,
                                   REML = TRUE)
          
          fit.stats <- r.squaredGLMM(
            fit.mod
          )
          
          LL <- logLik(
            fit.mod
          )[[1]]
          
          MLPE.list[[i + 1]] <-  fit.mod_REML
          
          cd.list[[i + 1]] <- as.matrix(cd)
          names(cd.list)[i + 1] <- "Distance"
          
          names(MLPE.list)[i + 1] <- "Distance"
          
          ROW <- nrow(gdist.inputs$ID)
          k <- 2
          
          k.list[[i + 1]] <- k
          names(k.list)[i + 1] <- 'Distance'
          
          
          
          n <- gdist.inputs$n.Pops
          aic <- (-2 * LL) + (2 * k)
          AICc <-
            # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
            (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
          
          if (GA.inputs$method == "AIC") {
            dist.obj <- -aic
          } else if (GA.inputs$method == "R2") {
            dist.obj <- fit.stats[[1]]
          } else {
            dist.obj <- LL[[1]]
          }
          
          Dist.AICc <- data.frame("Distance",
                                  dist.obj,
                                  k,
                                  aic,
                                  AICc,
                                  fit.stats[[1]],
                                  fit.stats[[2]],
                                  LL[[1]])
          colnames(Dist.AICc) <- c(
            "Surface",
            paste0("obj.func_", GA.inputs$method),
            'k',
            "AIC",
            "AICc",
            "R2m",
            "R2c",
            "LL"
          )
        }
        
        if (null_mod == TRUE) {
          dat <- gdist.inputs$df
          fit.null <- mlpe_rga(
            formula = .rga_null_model_formula(
              formula = gdist.inputs$formula,
              data = dat,
              fallback = gd ~ 1 + (1 | pop)
            ),
            data = dat,
            ZZ = gdist.inputs$ZZ,
            REML = FALSE
          )

          MLPE.list[['Null']] <- fit.null

          fit.stats <- r.squaredGLMM(fit.null)
          LL <- logLik(fit.null)
          ROW <- nrow(gdist.inputs$ID)
          k <- 1
          
          
          n <- gdist.inputs$n.Pops
          aic <- (-2 * LL) + (2 * k)
          AICc <-
            # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
            (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
          
          if (GA.inputs$method == "AIC") {
            null.obj <- -aic
          } else if (GA.inputs$method == "R2") {
            null.obj <- fit.stats[[1]]
          } else {
            null.obj <- LL[[1]]
          }
          
          Null.AICc <-
            data.frame("Null",
                       null.obj,
                       k,
                       aic,
                       AICc,
                       fit.stats[[1]],
                       fit.stats[[2]],
                       LL[[1]])
          colnames(Null.AICc) <-
            c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              'k',
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL"
            )
        }
      } # End Covariate
      
      # MLPE no Covariates ------------------------------------------------------
      if(is.null(gdist.inputs$covariates)) {
        
        # * Island GA ---------------------------------------------------------------
        
        
        if(isTRUE(GA.inputs$gaisl)) {
          
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              SCALE.vector(single.GA@solution, 1, GA.inputs$max.cat)
            
            # single.GA@solution <-
            #   single.GA@solution / min(single.GA@solution)
            df <- data.frame(id = terra::unique(r)[[1]], t(single.GA@solution))
            r <- terra::subst(r, from = df$id, to = df[[2]])
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,  "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = cd,
                genetic.dist = gdist.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "categorical",
                name = NAME,
                ID = gdist.inputs$ID,
                ZZ = gdist.inputs$ZZ
              )
            }
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # if (k.value == 1) {
            #   k <- 2
            # } else if (k.value == 2) {
            #   k <- GA.inputs$parm.type$n.parm[i] + 1
            # } else if (k.value == 3) {
            #   k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            # } else {
            #   k <- length(GA.inputs$layer.names[i]) + 1
            # }
            
            fit.mod <- mlpe_rga(formula = gdist.inputs$formula,
                                data = dat,
                                ZZ = gdist.inputs$ZZ,
                                REML = FALSE)
            
            fit.stats <- suppressWarnings(r.squaredGLMM(
              fit.mod
            ))
            
            LL <- logLik(
              fit.mod
            )[[1]]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + length(lme4::fixef(fit.mod)) - 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + length(lme4::fixef(fit.mod)) - 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + length(lme4::fixef(fit.mod)) - 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <-  MLPE.lmm2(
              resistance = cd,
              response = gdist.inputs$response,
              REML = TRUE,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # rm(single.GA, r)
            gc()
            
            # *-* Continuous -----------------------------------------------------------
            
          } else {
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method, "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = cd,
                genetic.dist = gdist.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "continuous",
                name = NAME,
                ID = gdist.inputs$ID,
                ZZ = gdist.inputs$ZZ
              )
              
              Plot.trans(
                PARM = single.GA@solution[-1],
                Resistance = GA.inputs$Resistance.stack[[i]],
                transformation = EQ,
                print.dir = GA.inputs$Plots.dir
              )
            }
            
            # fit.stats <-
            #   r.squaredGLMM(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            # 
            # 
            # aic <-
            #   AIC(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            # 
            # LL <-
            #   logLik(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            
            fit.mod <- mlpe_rga(formula = gdist.inputs$formula,
                                data = dat,
                                ZZ = gdist.inputs$ZZ,
                                REML = FALSE)
            
            fit.stats <- suppressWarnings(r.squaredGLMM(
              fit.mod
            ))
            
            LL <- logLik(
              fit.mod
            )[[1]]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + length(lme4::fixef(fit.mod)) - 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + length(lme4::fixef(fit.mod)) - 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + length(lme4::fixef(fit.mod)) - 1
            }
            
            MLPE.list[[i]] <-  MLPE.lmm2(
              resistance = cd,
              response = gdist.inputs$response,
              REML = TRUE,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # if (k.value == 1) {
            #   k <- 2
            # } else if (k.value == 2) {
            #   k <- GA.inputs$parm.type$n.parm[i] + 1
            # } else if (k.value == 3) {
            #   k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            # } else {
            #   k <- length(GA.inputs$layer.names[i]) + 1
            # }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
            # rm(single.GA, r)
            gc()
          } # Close gaisl cat-cont if else
        } else { # * Standard GA -------------------------------------------------------------
          
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              parallel = GA.inputs$parallel,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              SCALE.vector(single.GA@solution, 1, GA.inputs$max.cat)
            
            # single.GA@solution <-
            #   single.GA@solution / min(single.GA@solution)
            # df <- data.frame(id = terra::unique(r)[[1]], t(single.GA@solution))
            # r <- subs(r, df)
            df <- data.frame(id = terra::unique(GA.inputs$Resistance.stack[[i]])[[1]], t(single.GA@solution))
            r <- terra::subst(GA.inputs$Resistance.stack[[i]], from = df$id, to = df[[2]]) ## Modified 2019-03-26
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,  "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA, 
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = cd,
                genetic.dist = gdist.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "categorical",
                name = NAME,
                ID = gdist.inputs$ID,
                ZZ = gdist.inputs$ZZ
              )
            }
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # if (k.value == 1) {
            #   k <- 2
            # } else if (k.value == 2) {
            #   k <- GA.inputs$parm.type$n.parm[i] + 1
            # } else if (k.value == 3) {
            #   k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            # } else {
            #   k <- length(GA.inputs$layer.names[i]) + 1
            # }
            
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
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + length(lme4::fixef(fit.mod)) - 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + length(lme4::fixef(fit.mod)) - 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + length(lme4::fixef(fit.mod)) - 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # rm(single.GA, r)
            gc()
          } 
          else { # *-* Continuous ----------
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method, "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA, 
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = cd,
                genetic.dist = gdist.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "continuous",
                name = NAME,
                ID = gdist.inputs$ID,
                ZZ = gdist.inputs$ZZ
              )
              
              Plot.trans(
                PARM = single.GA@solution[-1],
                Resistance = GA.inputs$Resistance.stack[[i]],
                transformation = EQ,
                print.dir = GA.inputs$Plots.dir
              )
            }
            
            
            
            # fit.stats <-
            #   r.squaredGLMM(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            # 
            # 
            # aic <-
            #   AIC(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            # 
            # LL <-
            #   logLik(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            
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
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
            # rm(single.GA, r)
            gc()
          } # Close cat-cont if else
          
        } # Close GA vs. gaisl if-else
        
        if (dist_mod == TRUE) {
          r <- terra::classify(r, matrix(c(-Inf, Inf, 1), ncol = 3))
          names(r) <- "dist"
          cd <- Run_gdistance(gdist.inputs, r)
          
          dat <- gdist.inputs$df
          dat$cd <- scale(c(cd))
          
          write.table(
            as.matrix(cd),
            file = paste0(GA.inputs$Results.dir, 'Distance', "_", gdist.inputs$method, "_distMat.csv"),
            sep = ",",
            row.names = F,
            col.names = F
          )
          
          fit.mod <- mlpe_rga(formula = gdist.inputs$formula,
                              data = dat,
                              ZZ = gdist.inputs$ZZ,
                              REML = FALSE)
          
          fit.mod_REML <- mlpe_rga(formula = gdist.inputs$formula,
                                   data = dat,
                                   ZZ = gdist.inputs$ZZ,
                                   REML = TRUE)
          
          
          fit.stats <- r.squaredGLMM(
            fit.mod
          )
          
          LL <- logLik(
            fit.mod
          )[[1]]
          
          MLPE.list[[i + 1]] <-  fit.mod_REML
          
          cd.list[[i + 1]] <- as.matrix(cd)
          names(cd.list)[i + 1] <- "Distance"
          
          names(MLPE.list)[i + 1] <- "Distance"
          
          ROW <- nrow(gdist.inputs$ID)
          k <- 2
          
          k.list[[i + 1]] <- k
          names(k.list)[i + 1] <- 'Distance'
          
          n <- gdist.inputs$n.Pops
          aic <- (-2 * LL) + (2 * k)
          AICc <-
            # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
            (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
          
          if (GA.inputs$method == "AIC") {
            dist.obj <- -aic
          } else if (GA.inputs$method == "R2") {
            dist.obj <- fit.stats[[1]]
          } else {
            dist.obj <- LL[[1]]
          }
          
          Dist.AICc <- data.frame("Distance",
                                  dist.obj,
                                  k,
                                  aic,
                                  AICc,
                                  fit.stats[[1]],
                                  fit.stats[[2]],
                                  LL[[1]])
          colnames(Dist.AICc) <- c(
            "Surface",
            paste0("obj.func_", GA.inputs$method),
            'k',
            "AIC",
            "AICc",
            "R2m",
            "R2c",
            "LL"
          )
        }
        
        if (null_mod == TRUE) {
          dat <- gdist.inputs$df
          fit.null <- mlpe_rga(
            formula = .rga_null_model_formula(
              formula = gdist.inputs$formula,
              data = dat,
              fallback = gd ~ 1 + (1 | pop)
            ),
            data = dat,
            ZZ = gdist.inputs$ZZ,
            REML = FALSE
          )

          MLPE.list[['Null']] <- fit.null

          fit.stats <- r.squaredGLMM(fit.null)
          LL <- logLik(fit.null)
          ROW <- nrow(gdist.inputs$ID)
          k <- 1
          
          n <- gdist.inputs$n.Pops
          aic <- (-2 * LL) + (2 * k)
          AICc <-
            # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
            (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
          
          if (GA.inputs$method == "AIC") {
            null.obj <- -aic
          } else if (GA.inputs$method == "R2") {
            null.obj <- fit.stats[[1]]
          } else {
            null.obj <- LL[[1]]
          }
          
          Null.AICc <-
            data.frame("Null",
                       null.obj,
                       k,
                       aic,
                       AICc,
                       fit.stats[[1]],
                       fit.stats[[2]],
                       LL[[1]])
          colnames(Null.AICc) <-
            c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              'k',
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL"
            )
        }
        
      } # End no covariate
    } # End gdistance
    
    # >>> Julia <<< -------------------------------------------------
    if (!is.null(jl.inputs)) {
      
      # setwd(jl.inputs$JULIA_HOME)
      # MLPE with Covariates ----------------------------------------------------
      
      if(!is.null(jl.inputs$covariates)) { 
        
        # * Island GA ---------------------------------------------------------------
        
        
        if(isTRUE(GA.inputs$gaisl)) {
          stop("This feature is not currently supported")
          
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              SCALE.vector(single.GA@solution, 1, GA.inputs$max.cat)
            
            # single.GA@solution <-
            #   single.GA@solution / min(single.GA@solution)
            # df <- data.frame(id = terra::unique(r)[[1]], t(single.GA@solution))
            # r <- subs(r, df)
            df <- data.frame(id = terra::unique(GA.inputs$Resistance.stack[[i]])[[1]], t(single.GA@solution))
            r <- terra::subst(GA.inputs$Resistance.stack[[i]], from = df$id, to = df[[2]]) ## Modified 2019-03-26
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,  "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = cd,
                genetic.dist = gdist.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "categorical",
                name = NAME,
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
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 
                length(lme4::fixef(fit.mod)) - 1
              
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + 
                length(GA.inputs$layer.names) + 
                length(lme4::fixef(fit.mod)) - 1
              
            } else {
              k <- length(GA.inputs$layer.names[i]) + 
                length(lme4::fixef(fit.mod)) - 1
              
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <- fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # rm(single.GA, r)
            gc()
            
            # *-* Continuous -----------------------------------------------------------
            
          } else {
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method, "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = cd,
                genetic.dist = gdist.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "continuous",
                name = NAME,
                ID = gdist.inputs$ID,
                ZZ = gdist.inputs$ZZ
              )
              
              Plot.trans(
                PARM = single.GA@solution[-1],
                Resistance = GA.inputs$Resistance.stack[[i]],
                transformation = EQ,
                print.dir = GA.inputs$Plots.dir
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
            
            fit.stats <- r.squaredGLMM(
              fit.mod
            )
            
            LL <- logLik(
              fit.mod
            )[[1]]
            
            MLPE.list[[i + 1]] <- fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 
                length(lme4::fixef(fit.mod)) - 1
              
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + 
                length(GA.inputs$layer.names) + 
                length(lme4::fixef(fit.mod)) - 1
              
            } else {
              k <- length(GA.inputs$layer.names[i]) + 
                length(lme4::fixef(fit.mod)) - 1
              
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
            # rm(single.GA, r)
            gc()
          } # Close gaisl cat-cont if else
        } else { # * Standard GA -------------------------------------------------------------
          
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              jl.inputs = jl.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              parallel = GA.inputs$parallel,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            ## Turned off 11/14/2023
            # single.GA@solution <-
            #   SCALE.vector(single.GA@solution, 1, GA.inputs$max.cat)
            
            ##
            ga.p <- single.GA@solution
            parm <- ga.p / min(ga.p)
            if(max(parm) > GA.inputs$max.cat){
              parm <- SCALE.vector(parm, 1, GA.inputs$max.cat)
            }
            
            single.GA@solution <- parm
            ##
            
            df <- data.frame(id = terra::unique(GA.inputs$Resistance.stack[[i]])[[1]], t(single.GA@solution))
            r <- terra::subst(GA.inputs$Resistance.stack[[i]], from = df$id, to = df[[2]]) ## Modified 2019-03-26
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA, 
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = dat$cd,
                genetic.dist = jl.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "categorical",
                name = NAME,
                ID = jl.inputs$ID,
                ZZ = jl.inputs$ZZ
              )
            }
            
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
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 
                length(lme4::fixef(fit.mod)) - 1
              
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + 
                length(GA.inputs$layer.names) + 
                length(lme4::fixef(fit.mod)) - 1
              
            } else {
              k <- length(GA.inputs$layer.names[i]) + 
                length(lme4::fixef(fit.mod)) - 1
              
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- jl.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <- fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # rm(single.GA, r)
            gc()
          } 
          else { # *-* Continuous ----------
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              jl.inputs = jl.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA, 
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = dat$cd,
                genetic.dist = jl.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "continuous",
                name = NAME,
                ID = jl.inputs$ID,
                ZZ = jl.inputs$ZZ
              )
              
              Plot.trans(
                PARM = single.GA@solution[-1],
                Resistance = GA.inputs$Resistance.stack[[i]],
                transformation = EQ,
                print.dir = GA.inputs$Plots.dir
              )
            }
            
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
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + length(lme4::fixef(fit.mod)) - 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + length(lme4::fixef(fit.mod)) - 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + length(lme4::fixef(fit.mod)) - 1
            }
            
            MLPE.list[[i]] <- fit.mod_REML
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- jl.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
            # rm(single.GA, r)
            gc()
          } # Close cat-cont if else
          
        } # Close GA vs. gaisl if-else
        
        if (dist_mod == TRUE) {
          r <- terra::classify(r, matrix(c(-Inf, Inf, 1), ncol = 3))
          names(r) <- "dist"
          cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
          cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
          dat <- jl.inputs$df
          dat$cd <- cd.l
          
          write.table(
            cd,
            file = paste0(GA.inputs$Results.dir, "Distance", "_jlResistMat.csv"),
            sep = ",",
            row.names = F,
            col.names = F
          )
          
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
          
          MLPE.list[[i + 1]] <-  fit.mod_REML
          
          cd.list[[i + 1]] <- cd
          names(cd.list)[i + 1] <- "Distance"
          
          names(MLPE.list)[i + 1] <- "Distance"
          
          ROW <- nrow(jl.inputs$ID)
          k <- length(lme4::fixef(fit.mod))
          
          k.list[[i + 1]] <- k
          names(k.list)[i + 1] <- 'Distance'
          
          
          n <- jl.inputs$n.Pops
          aic <- (-2 * LL) + (2 * k)
          AICc <-
            # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
            (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
          
          if (GA.inputs$method == "AIC") {
            dist.obj <- -aic
          } else if (GA.inputs$method == "R2") {
            dist.obj <- fit.stats[[1]]
          } else {
            dist.obj <- LL[[1]]
          }
          
          Dist.AICc <- data.frame("Distance",
                                  dist.obj,
                                  k,
                                  aic,
                                  AICc,
                                  fit.stats[[1]],
                                  fit.stats[[2]],
                                  LL[[1]])
          colnames(Dist.AICc) <- c(
            "Surface",
            paste0("obj.func_", GA.inputs$method),
            'k',
            "AIC",
            "AICc",
            "R2m",
            "R2c",
            "LL"
          )
        }
        
        if (null_mod == TRUE) {
          dat <- jl.inputs$df
          fit.null <- mlpe_rga(
            formula = .rga_null_model_formula(
              formula = jl.inputs$formula,
              data = dat,
              fallback = gd ~ 1 + (1 | pop)
            ),
            data = dat,
            ZZ = jl.inputs$ZZ,
            REML = FALSE
          )

          MLPE.list[['Null']] <- fit.null

          fit.stats <- r.squaredGLMM(fit.null)
          LL <- logLik(fit.null)
          ROW <- nrow(jl.inputs$ID)
          k <- 1
          
          
          n <- jl.inputs$n.Pops
          aic <- (-2 * LL) + (2 * k)
          AICc <-
            # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
            (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
          
          if (GA.inputs$method == "AIC") {
            null.obj <- -aic
          } else if (GA.inputs$method == "R2") {
            null.obj <- fit.stats[[1]]
          } else {
            null.obj <- LL[[1]]
          }
          
          Null.AICc <-
            data.frame("Null",
                       null.obj,
                       k,
                       aic,
                       AICc,
                       fit.stats[[1]],
                       fit.stats[[2]],
                       LL[[1]])
          colnames(Null.AICc) <-
            c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              'k',
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL"
            )
        }
      } # End Covariate
      
      
      # MLPE no Covariates ------------------------------------------------------
      if(is.null(jl.inputs$covariates)) { 
        
        
        # * Island GA ---------------------------------------------------------------
        
        if(isTRUE(GA.inputs$gaisl)) {
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              jl.inputs = jl.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              SCALE.vector(single.GA@solution, 1, GA.inputs$max.cat)
            
            # single.GA@solution <-
            #   single.GA@solution / min(single.GA@solution)
            # df <- data.frame(id = terra::unique(r)[[1]], t(single.GA@solution))
            # r <- subs(r, df)
            df <- data.frame(id = terra::unique(GA.inputs$Resistance.stack[[i]])[[1]], t(single.GA@solution))
            r <- terra::subst(GA.inputs$Resistance.stack[[i]], from = df$id, to = df[[2]]) ## Modified 2019-03-26
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = dat$cd,
                genetic.dist = jl.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "categorical",
                name = NAME,
                ID = jl.inputs$ID,
                ZZ = jl.inputs$ZZ
              )
            }
            
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
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = dat$cd,
            #     response = jl.inputs$response,
            #     REML = F,
            #     ID = jl.inputs$ID,
            #     ZZ = jl.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = dat$cd,
            #     response = jl.inputs$response,
            #     REML = F,
            #     ID = jl.inputs$ID,
            #     ZZ = jl.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = dat$cd,
            #     response = jl.inputs$response,
            #     REML = F,
            #     ID = jl.inputs$ID,
            #     ZZ = jl.inputs$ZZ
            #   )
            # )
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- jl.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- cd
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
          } else { # *-* Continuous ---------------
            
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              jl.inputs = jl.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = dat$cd,
                genetic.dist = jl.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "continuous",
                name = NAME,
                ID = jl.inputs$ID,
                ZZ = jl.inputs$ZZ
              )
              
              Plot.trans(
                PARM = single.GA@solution[-1],
                Resistance = GA.inputs$Resistance.stack[[i]],
                transformation = EQ,
                print.dir = GA.inputs$Plots.dir
              )
            }
            
            
            
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
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- cd
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- jl.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
          } # Close if-else
          
          if (dist_mod == TRUE) {
            r <- terra::classify(r, matrix(c(-Inf, Inf, 1), ncol = 3))
            names(r) <- "dist"
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, "Distance", "_jlResistMat.csv"),
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            fit.mod <- mlpe_rga(formula = jl.inputs$formula,
                                data = dat,
                                ZZ = jl.inputs$ZZ,
                                REML = FALSE)
            fit.mod_REML <- mlpe_rga(formula = jl.inputs$formula,
                                     data = dat,
                                     ZZ = jl.inputs$ZZ,
                                     REML = TRUE)
            
            
            fit.stats <- r.squaredGLMM(
              fit.mod
            )
            
            LL <- logLik(
              fit.mod
            )[[1]]
            
            MLPE.list[[i + 1]] <-  fit.mod_REML
            
            cd.list[[i + 1]] <- cd
            names(cd.list)[i + 1] <- "Distance"
            
            names(MLPE.list)[i + 1] <- "Distance"
            
            ROW <- nrow(jl.inputs$ID)
            k <- 2
            
            k.list[[i + 1]] <- k
            names(k.list)[i + 1] <- 'Distance'
            
            
            n <- jl.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            if (GA.inputs$method == "AIC") {
              dist.obj <- -aic
            } else if (GA.inputs$method == "R2") {
              dist.obj <- fit.stats[[1]]
            } else {
              dist.obj <- LL[[1]]
            }
            
            Dist.AICc <- data.frame("Distance",
                                    dist.obj,
                                    k,
                                    aic,
                                    AICc,
                                    fit.stats[[1]],
                                    fit.stats[[2]],
                                    LL[[1]])
            colnames(Dist.AICc) <- c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              'k',
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL"
            )
          }
          
          if (null_mod == TRUE) {
            dat <- jl.inputs$df
            fit.null <- mlpe_rga(
              formula = .rga_null_model_formula(
                formula = jl.inputs$formula,
                data = dat,
                fallback = gd ~ 1 + (1 | pop)
              ),
              data = dat,
              ZZ = jl.inputs$ZZ,
              REML = FALSE
            )

            MLPE.list[['Null']] <- fit.null

          fit.stats <- r.squaredGLMM(fit.null)
            LL <- logLik(fit.null)
            ROW <- nrow(jl.inputs$ID)
            k <- 1
            
            n <- jl.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            if (GA.inputs$method == "AIC") {
              null.obj <- -aic
            } else if (GA.inputs$method == "R2") {
              null.obj <- fit.stats[[1]]
            } else {
              null.obj <- LL[[1]]
            }
            
            Null.AICc <-
              data.frame("Null",
                         null.obj,
                         k,
                         aic,
                         AICc,
                         fit.stats[[1]],
                         fit.stats[[2]],
                         LL[[1]])
            colnames(Null.AICc) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL"
              )
          }
          
          
        } else { # * Standard GA -----------------
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              jl.inputs = jl.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            ## Turned off 11/14/2023
            # single.GA@solution <-
            #   SCALE.vector(single.GA@solution, 1, GA.inputs$max.cat)
            
            ##
            ga.p <- single.GA@solution
            parm <- ga.p / min(ga.p)
            if(max(parm) > GA.inputs$max.cat){
              parm <- SCALE.vector(parm, 1, GA.inputs$max.cat)
            }
            
            single.GA@solution <- parm
            ##
            
            df <- data.frame(id = terra::unique(GA.inputs$Resistance.stack[[i]])[[1]], t(single.GA@solution))
            r <- terra::subst(GA.inputs$Resistance.stack[[i]], from = df$id, to = df[[2]]) ## Modified 2019-03-26
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = dat$cd,
                genetic.dist = jl.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "categorical",
                name = NAME,
                ID = jl.inputs$ID,
                ZZ = jl.inputs$ZZ
              )
            }
            
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
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = dat$cd,
            #     response = jl.inputs$response,
            #     REML = F,
            #     ID = jl.inputs$ID,
            #     ZZ = jl.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = dat$cd,
            #     response = jl.inputs$response,
            #     REML = F,
            #     ID = jl.inputs$ID,
            #     ZZ = jl.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = dat$cd,
            #     response = jl.inputs$response,
            #     REML = F,
            #     ID = jl.inputs$ID,
            #     ZZ = jl.inputs$ZZ
            #   )
            # )
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- jl.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- cd
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
          } else { # *-* Continuous ---------------
            
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = .rga_wrap_raster_for_parallel(r, GA.inputs$parallel),
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.worker.inputs,
              jl.inputs = jl.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            if(isTRUE(diagnostic_plots)){
              Diagnostic.Plots(
                resistance.mat = dat$cd,
                genetic.dist = jl.inputs$response,
                plot.dir = GA.inputs$Plots.dir,
                type = "continuous",
                name = NAME,
                ID = jl.inputs$ID,
                ZZ = jl.inputs$ZZ
              )
              
              Plot.trans(
                PARM = single.GA@solution[-1],
                Resistance = GA.inputs$Resistance.stack[[i]],
                transformation = EQ,
                print.dir = GA.inputs$Plots.dir
              )
            }
            
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
            
            # fit.stats <-
            #   r.squaredGLMM(
            #     MLPE.lmm(
            #       resistance = dat$cd,
            #       pairwise.genetic = jl.inputs$response,
            #       REML = F,
            #       ID = jl.inputs$ID,
            #       ZZ = jl.inputs$ZZ
            #     )
            #   )
            # 
            # 
            # aic <-
            #   AIC(
            #     MLPE.lmm(
            #       resistance = dat$cd,
            #       pairwise.genetic = jl.inputs$response,
            #       REML = F,
            #       ID = jl.inputs$ID,
            #       ZZ = jl.inputs$ZZ
            #     )
            #   )
            # 
            # LL <-
            #   logLik(
            #     MLPE.lmm(
            #       resistance = dat$cd,
            #       pairwise.genetic = jl.inputs$response,
            #       REML = F,
            #       ID = jl.inputs$ID,
            #       ZZ = jl.inputs$ZZ
            #     )
            #   )
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- cd
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- jl.inputs$n.Pops
            
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
          } # Close  cat-cont ifelse
          
          if (dist_mod == TRUE) {
            r <- terra::classify(r, matrix(c(-Inf, Inf, 1), ncol = 3))
            names(r) <- "dist"
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, "Distance", "_jlResistMat.csv"),
              sep = ",",
              row.names = F,
              col.names = F
            )
            
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
            
            MLPE.list[[i + 1]] <-  fit.mod_REML
            
            cd.list[[i + 1]] <- cd
            names(cd.list)[i + 1] <- "Distance"
            
            names(MLPE.list)[i + 1] <- "Distance"
            
            ROW <- nrow(jl.inputs$ID)
            k <- 2
            
            k.list[[i + 1]] <- k
            names(k.list)[i + 1] <- 'Distance'
            
            n <- jl.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            if (GA.inputs$method == "AIC") {
              dist.obj <- -aic
            } else if (GA.inputs$method == "R2") {
              dist.obj <- fit.stats[[1]]
            } else {
              dist.obj <- LL[[1]]
            }
            
            Dist.AICc <- data.frame("Distance",
                                    dist.obj,
                                    k,
                                    aic,
                                    AICc,
                                    fit.stats[[1]],
                                    fit.stats[[2]],
                                    LL[[1]])
            colnames(Dist.AICc) <- c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              'k',
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL"
            )
          }
          
          if (null_mod == TRUE) {
            dat <- jl.inputs$df
            fit.null <- mlpe_rga(
              formula = .rga_null_model_formula(
                formula = jl.inputs$formula,
                data = dat,
                fallback = gd ~ 1 + (1 | pop)
              ),
              data = dat,
              ZZ = jl.inputs$ZZ,
              REML = FALSE
            )

            MLPE.list[['Null']] <- fit.null

          fit.stats <- r.squaredGLMM(fit.null)
            LL <- logLik(fit.null)
            ROW <- nrow(jl.inputs$ID)
            k <- 1
            
            
            n <- jl.inputs$n.Pops
            aic <- (-2 * LL) + (2 * k)
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            if (GA.inputs$method == "AIC") {
              null.obj <- -aic
            } else if (GA.inputs$method == "R2") {
              null.obj <- fit.stats[[1]]
            } else {
              null.obj <- LL[[1]]
            }
            
            Null.AICc <-
              data.frame("Null",
                         null.obj,
                         k,
                         aic,
                         AICc,
                         fit.stats[[1]],
                         fit.stats[[2]],
                         LL[[1]])
            colnames(Null.AICc) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL"
              )
          }
          
        } # End Island-Standard ifelse
      } # End no covariates
    } # End julia
  } # Close raster loop
  
  
  # Optimization summary ----------------------------------------------------
  
  # Make results data frame
  Results.cat <- data.frame()
  Results.cont <- data.frame()
  # cnt1<-0
  # cnt2<-0
  for (i in 1:GA.inputs$n.layers) {
    if (GA.inputs$surface.type[i] == 'cat') {
      #     cnt1 <- cnt1+1
      #     RS <- data.frame(GA.inputs$layer.names[i], -(RESULTS.cat[[i]]@fitnessValue),RESULTS[[i]]@solution)
      Results.cat <- do.call(rbind.fill, RESULTS.cat)
    } else {
      #   cnt2 <-cnt2+1
      #   RS <- data.frame(GA.inputs$layer.names[i], -(RESULTS.cont[[i]]@fitnessValue), Cont.Param(RESULTS[[i]]@solution))
      Results.cont <- do.call(rbind, RESULTS.cont)
    }
  }
  
  
  n.pops <- if (!is.null(gdist.inputs)) gdist.inputs$n.Pops else jl.inputs$n.Pops

  # Compile results into tables
  cat("\n")
  cat("\n")
  if (nrow(Results.cat) > 0) {
    Features <- array()
    for (i in 1:ncol(Results.cat) - 8) {
      feature <- paste0("Feature", i)
      Features[i] <- feature
    }
    colnames(Results.cat) <-
      c(
        "Surface",
        paste0("obj.func_", GA.inputs$method),
        'k',
        "AIC",
        "AICc",
        "R2m",
        "R2c",
        "LL",
        Features
      )
    Results.cat <- .rga_update_result_table_models(
      result_table = Results.cat,
      fit_list = MLPE.list,
      GA.inputs = GA.inputs,
      n_pops = n.pops
    )
    Results.cat <-  Results.cat[order(Results.cat$AICc), ]
    write.table(
      Results.cat,
      paste0(GA.inputs$Results.dir, "CategoricalResults.csv"),
      sep = ",",
      col.names = T,
      row.names = F
    )
  }
  
  if (ncol(Results.cont) > 0) {
    colnames(Results.cont) <-
      c(
        "Surface",
        paste0("obj.func_", GA.inputs$method),
        'k',
        "AIC",
        "AICc",
        "R2m",
        "R2c",
        "LL",
        "Equation",
        "shape",
        "max"
      )
    Results.cont <- .rga_update_result_table_models(
      result_table = Results.cont,
      fit_list = MLPE.list,
      GA.inputs = GA.inputs,
      n_pops = n.pops
    )
    Results.cont <- Results.cont[order(Results.cont$AICc), ]
    write.table(
      Results.cont,
      paste0(GA.inputs$Results.dir, "ContinuousResults.csv"),
      sep = ",",
      col.names = T,
      row.names = F
    )
  }
  
  # Full Results
  summary.cols <- c(
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

  if (nrow(Results.cat) > 0 & nrow(Results.cont) > 0) {
    Results.All <- rbind(
      Results.cat[, summary.cols, drop = FALSE],
      Results.cont[, summary.cols, drop = FALSE]
    )
  } else if (nrow(Results.cat) < 1 & nrow(Results.cont) > 0) {
    Results.All <- Results.cont[, summary.cols, drop = FALSE]
  } else {
    Results.All <- Results.cat[, summary.cols, drop = FALSE]
  }
  
  if (dist_mod == TRUE) {
    Dist.AICc$Fixed.Effects <- NA_character_
    Dist.AICc <- Dist.AICc[, summary.cols, drop = FALSE]
    Results.All <- rbind(Results.All, Dist.AICc)
  }

  if (null_mod == TRUE) {
    Null.AICc$Fixed.Effects <- NA_character_
    Null.AICc <- Null.AICc[, summary.cols, drop = FALSE]
    Results.All <- rbind(Results.All, Null.AICc)
  }
  
  Results.All <- .rga_update_result_table_models(
    result_table = Results.All,
    fit_list = MLPE.list,
    GA.inputs = GA.inputs,
    n_pops = n.pops
  )
  Results.All <- Results.All[order(Results.All$AICc), ]
  
  cat("\n")
  cat("\n")
  write.table(
    Results.All,
    paste0(GA.inputs$Results.dir, "All_Results_Table.csv"),
    
    sep = ",",
    col.names = T,
    row.names = F
  )
  
  # if(!is.null(gdist.inputs$covariates)) { 
  #   MLPE.results <- NULL
  # } else {
  # Get parameter estimates
  if (!is.null(jl.inputs)) {
    MLPE.results <- MLPE.lmm_coef(
      formula = jl.inputs$formula,
      inputs = jl.inputs$df,
      resistance = GA.inputs$Results.dir,
      genetic.dist = jl.inputs$response,
      out.dir = GA.inputs$Results.dir,
      method = "jl",
      ID = jl.inputs$ID,
      ZZ = jl.inputs$ZZ
    )
  } else {
    MLPE.results <- MLPE.lmm_coef(
      resistance = GA.inputs$Results.dir,
      genetic.dist = gdist.inputs$response,
      out.dir = GA.inputs$Results.dir,
      method = "gd",
      ID = gdist.inputs$ID,
      ZZ = gdist.inputs$ZZ
    )
  } 
  # } ## End covariate if-else
  
  k.list <- Results.All[, c("Surface", "k", "Fixed.Effects"), drop = FALSE]
  colnames(k.list) <- c("surface", "k", "fixed.effects")
  
  rt <- proc.time()[3] - t1
  # Full Results
  if (nrow(Results.cat) > 0 & nrow(Results.cont) > 0) {
    RESULTS <-
      list(
        ContinuousResults = Results.cont,
        CategoricalResults = Results.cat,
        AICc = Results.All,
        MLPE = MLPE.results,
        Run.Time = rt,
        MLPE.list = MLPE.list,
        cd = cd.list,
        k = k.list,
        ga = ga.list
      )
    
  } else if (nrow(Results.cat) < 1 & nrow(Results.cont) > 0) {
    RESULTS <-
      list(
        ContinuousResults = Results.cont,
        CategoricalResults = NULL,
        AICc = Results.All,
        MLPE = MLPE.results,
        Run.Time = rt,
        MLPE.list = MLPE.list,
        cd = cd.list,
        k = k.list,
        ga = ga.list
      )
    
  } else if (nrow(Results.cat) > 0 & nrow(Results.cont) < 1) {
    RESULTS <-
      list(
        ContinuousResults = NULL,
        CategoricalResults = Results.cat,
        AICc = Results.All,
        MLPE = MLPE.results,
        Run.Time = rt,
        MLPE.list = MLPE.list,
        cd = cd.list,
        k = k.list,
        ga = ga.list
      )
  } else {
    RESULTS <-
      list(
        ContinuousResults = NULL,
        CategoricalResults = NULL,
        AICc = Results.All,
        MLPE = MLPE.results,
        Run.Time = rt,
        MLPE.list = MLPE.list,
        cd = cd.list,
        k = k.list,
        ga = ga.list
      )
  }
  # rm(single.GA, r)
  setwd(wd)
  gc()
  RESULTS <- resga_add_class(RESULTS, "resga_ss_optim")
  return(RESULTS)
}



