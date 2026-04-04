#' Analyze all combinations
#'
#' A wrapper function to run both \code{\link[ResistanceGA2]{SS_optim}} and
#' \code{\link[ResistanceGA2]{MS_optim}} to optimize all combinations of
#' resistance surfaces with the Genetic Algorithms package \pkg{GA}. Following
#' optimization, \code{\link[ResistanceGA2]{Resist.boot}} is run to conduct a
#' bootstrap analysis.
#'
#' @param gdist.inputs Object created from running
#'   \code{\link[ResistanceGA2]{gdist.prep}}. Supply when optimizing with
#'   \pkg{gdistance}.
#' @param jl.inputs Object created from running
#'   \code{\link[ResistanceGA2]{jl.prep}}. Supply when optimizing with
#'   Circuitscape in Julia.
#' @param GA.inputs Object created from running
#'   \code{\link[ResistanceGA2]{GA.prep}}. Be sure that
#'   \code{Results.dir = "all.comb"} was used when creating this object.
#' @param results.dir Directory to write and save analysis results. If the
#'   directory does not exist, it will be created automatically. Any existing
#'   files located in this directory may be deleted after confirmation.
#' @param max.combination The maximum number of surfaces to include in the all combinations analysis (Default = 4). Alternatively, specify a vector with the minimum and maximum number of surfaces to combine (e.g., c(2,4). If the minimum > 1, then the single surface optimization will be skipped.
#' @param iters Number of bootstrap iterations to be conducted (Default = 1000)
#' @param sample.prop Proportion of observations to be sampled each iteration (Default = 0.75)
#' @param replicate The number of times to replicate the GA optimization process for each surface (Default = 1)
#' @param nlm (NOT CURRENTLY IMPLEMENTED) Logical, if TRUE, the final step of optimization will use nlm to fine-tune parameter estimates. This may lead to overfitting in some cases. (Default = FALSE)
#' @param dist_mod Logical, if TRUE, a Distance model will be calculated and added to the output table (default = TRUE)
#' @param null_mod Logical, if TRUE, an intercept-only model will be calculated and added to the output table (Default = TRUE)
#' @param ... Additional arguments (Not currently used)

#' @return An object of class \code{resga_all_comb}. This function optimizes
#'   resistance surfaces in isolation using
#'   \code{\link[ResistanceGA2]{SS_optim}}, followed by multisurface
#'   optimization using \code{\link[ResistanceGA2]{MS_optim}}, and then
#'   conducts a bootstrap analysis.
#' @usage all_comb(gdist.inputs = NULL,
#'                 jl.inputs = NULL,
#'                 GA.inputs,
#'                 results.dir,
#'                 max.combination = 4,
#'                 iters = 1000,
#'                 replicate = 1,
#'                 sample.prop = 0.75,
#'                 nlm = FALSE,
#'                 dist_mod = TRUE,
#'                 null_mod = TRUE,
#'                 ...)
#' @author Bill Peterman <Peterman.73@@osu.edu>
#' @export
#' 
#' @examples
#' \dontrun{
#' pts <- terra::vect(sample_pops[[1]], type = "points")
#' gd <- lower(Dc_list[[1]])
#'
#' gdist.inputs <- gdist.prep(
#'   n.Pops = nrow(sample_pops[[1]]),
#'   response = gd,
#'   samples = pts
#' )
#'
#' GA.inputs <- GA.prep(
#'   raster = raster_orig,
#'   Results.dir = "all.comb",
#'   quiet = TRUE
#' )
#'
#' results_dir <- file.path(tempdir(), "ResistanceGA2-all-comb")
#'
#' combo.out <- all_comb(
#'   gdist.inputs = gdist.inputs,
#'   GA.inputs = GA.inputs,
#'   results.dir = results_dir,
#'   max.combination = 2,
#'   iters = 10,
#'   replicate = 1
#' )
#' }

all_comb <- function(gdist.inputs = NULL,
                     jl.inputs = NULL,
                     GA.inputs,
                     results.dir,
                     max.combination = 4,
                     iters = 1000,
                     replicate = 1,
                     sample.prop = 0.75,
                     nlm = FALSE,
                     dist_mod = TRUE,
                     null_mod = TRUE,
                     ...) {
  
  if(!exists('results.dir')) 
    return(cat("ERROR: An empty results directory must be specified"))
  
  if(!exists('gdist.inputs')) 
    return(cat("ERROR: Please specify gdist.inputs"))
  
  # if(!is.null('CS.inputs')) 
  #   return(cat("ERROR: All combinations analysis cannot currently be done using CIRCUITSCAPE"))
  
  if(!exists('GA.inputs')) 
    return(cat("ERROR: Please specify GA.inputs"))
  
  if(!is.null(GA.inputs$Results.dir) & 
     !is.null(GA.inputs$Write.dir) &
     !is.null(GA.inputs$Plots.dir)) {
    return(cat("ERROR: Please correctly specify the `Results.dir` as 'all.comb' when running GA.prep"))
  }
  
  if(length(max.combination) > 2) {
    return(cat("ERROR: Please specify either a single value or a vector with the minimum and maximum value"))
  }
  
  results.dir <- paste0(
    sub("[/\\\\]+$", "", normalizePath(results.dir, winslash = "/", mustWork = FALSE)),
    "/"
  )
  
  if(!dir.exists(results.dir)) {
    dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)
    if(!dir.exists(results.dir)) {
      stop("Failed to create 'results.dir': ", results.dir)
    }
    message("Created 'results.dir': ", results.dir)
  }
  
  dir.files <- list.files(results.dir)
  
  dotargs <- list(...)
  if(is.null(dotargs$cluster)) dotargs$cluster <- FALSE
  if(parallel::detectCores() > 16) dotargs$cluster <- TRUE
  
  if(isTRUE(dotargs$cluster) && length(dir.files) != 0) {
    unlink(dir(results.dir, 
               full.names = TRUE),
           recursive = TRUE,
           force = T)
  } 
  
  dir.files <- list.files(results.dir)
  
  if(length(dir.files) != 0) {
    q <- yn.question(cat(
      paste0("This function is about to delete all files and folders in '", results.dir, "'"),
      '\n', '\n',
      paste0("Do you want to proceed? Select 1 (Yes), 2 (No), or 3 (create a new subdirectory),  then press [Enter]")))
    
    # if(q == FALSE) return(cat("Function stopped"))
    
    if(is.na(q)) { # Create subdir
      dir.NAME <- floor(as.numeric(as.POSIXct(Sys.time())))
      dir.create(path = paste0(results.dir, "all_comb_", dir.NAME))
      results.dir <- paste0(results.dir, "all_comb_", dir.NAME, "/")
    } else if(q == FALSE) { # Stop function
      return(cat("Function stopped"))
    } else { # Remove exisiting folder
      unlink(dir(results.dir, 
                 full.names = TRUE),
             recursive = TRUE,
             force = T
      )
    }
  }
  
  
  
  # if(length(max.combination) == 2) {
  #   if(max.combination[2] > GA.inputs$n.layers) {
  #     return(cat("ERROR: Please specify a maximum combination that is less than or equal to the number of raster layers in the analysis"))
  #   }
  # }
  
  
  # Create combination list -------------------------------------------------
## Check on creation of combination list!!!
  ## Errored when trying to run select combinations.
    mc <- max.combination
  
  if(length(max.combination) == 2) {
    if(mc[1] == 1) {
      min.combination <- 2
      max.combination <- mc[2]
      ss <- TRUE
    } else {
      min.combination <- mc[1]
      max.combination <- mc[2]
      ss <- FALSE
    } 
  } else {
    min.combination <- 2
    ss <- TRUE
  }
  
  if(max.combination > GA.inputs$n.layers) {
    max.combination <- GA.inputs$n.layers
  }
  
  comb.list <- vector(mode = "list", length = (max.combination - 1))
  
  
  list.count <- 0
  surface.count <- 0
  for(i in min.combination:max.combination) {
    list.count <- list.count + 1
    comb.list[[list.count]] <- t(combn(1:GA.inputs$n.layers, i))
    if(is.null(nrow(comb.list[[list.count]]))) {
      n.comb <- 1
    } else {
      n.comb <- nrow(comb.list[[list.count]])
    }
    surface.count <- surface.count + n.comb
  }
  
  all.combs <- list()
  comb.names <- list()
  row.index <- 0
  for(i in 1:length(comb.list)){
    combs <- comb.list[[i]]
    
    if(is.null(nrow(comb.list[[i]]))) {
      t.combs <- 1
    } else {
      t.combs <- nrow(comb.list[[i]])
    }
    
    for(j in 1:t.combs) {
      row.index <- row.index + 1
      all.combs[[row.index]] <- combs[j,]
      c.names <- GA.inputs$layer.names[combs[j,]]
      comb.names[[row.index]] <- paste(c.names, collapse = ".")
    }
  }
  
  GA.input_orig <- GA.inputs
  
  Results <- vector(mode = 'list', length = replicate)
  # Begin Replicate Loop --------------------------------------------------
  for(i in 1:replicate){
    # Skip if min combination > 1
    if(ss == FALSE) {
      ss.results <- NULL
      AICc.tab <- NULL
      dir.create(paste0(results.dir,'rep_',i))
    } else {  # Do single surface optimization
      dir.create(paste0(results.dir,'rep_',i))
      dir.create(paste0(results.dir,'rep_',i, "/", "single.surface"))
      dir.create(paste0(results.dir,'rep_',i, "/", "single.surface/", "Plots"))
      
      # Single-surface optimization -------------------------------------------
      GA.inputs$Plots.dir <- paste0(results.dir,
                                    'rep_',i, 
                                    "/",
                                    "single.surface/",
                                    "Plots/")
      
      GA.inputs$Results.dir <- paste0(results.dir,
                                      'rep_',i, 
                                      "/", 
                                      "single.surface/")
      
      if(!is.null(gdist.inputs)) {
        ss.results <- SS_optim(gdist.inputs = gdist.inputs,
                               GA.inputs = GA.inputs,
                               dist_mod = dist_mod,
                               null_mod = null_mod)
      } else if(!is.null(jl.inputs)) {
        ss.results <- SS_optim(jl.inputs = jl.inputs,
                               GA.inputs = GA.inputs,
                               dist_mod = dist_mod,
                               null_mod = null_mod)
      } else {
        stop("`gdist.inputs` or `jl.inputs` must be specified!!!")
      }
      
      AICc.tab <- ss.results$AICc
    }
    
    if(is.null(ss.results)) {
      best.list <- list()
    } else if(isTRUE(GA.input_orig$gaisl)) {
      best.list <- lapply(ss.results$ga, function(x) x@solutions) # Extract best individuals
      best.list <- lapply(best.list, function(x) plyr::ldply(x))
    } else {
      best.list <- lapply(ss.results$ga, function(x) x@population) # Extract best individuals
    }
    
    
    # Multisurface optimization -----------------------------------------------
    
    ms.cd <- vector(mode = 'list',
                    length = length(all.combs))
    
    ms.k <- vector(mode = 'list',
                   length = length(all.combs))
    
    AICc.tab_list <- vector(mode = 'list',
                            length = length(all.combs))
    
    ms.results <- vector(mode = "list", length = length(all.combs))
    
    if(is.null(ss.results)) {
      n_ss.cd <- 0
      all.cd <- ms.cd
    } else {
      n_ss.cd <- length(ss.results$cd)
      all.cd <- c(ss.results$cd, ms.cd)
    }
    
    for(j in 1:length(all.combs)) {
        dir.create(paste0(results.dir,'rep_',i, "/", comb.names[[j]]))
        # dir.create(paste0(results.dir,'rep_',i, "/", comb.names[[j]],"/Plots"))
        
        # Select raster surfaces
        raster.comb <- terra::subset(GA.input_orig$Resistance.stack, all.combs[[j]])

        if(!is.null(GA.input_orig$inputs$select.trans)) {
          s.trans <- GA.input_orig$inputs$select.trans[all.combs[[j]]]
        } else {
          s.trans <- GA.input_orig$inputs$select.trans
        }
        
        ## Make suggestions
        suggest.sample <- sample(GA.input_orig$pop.size, floor(GA.input_orig$pop.size * 0.5), replace = F)
        suggest <- vector(mode = 'list', length = length(all.combs[[j]]))
        
        if(length(best.list) != 0){
          for(s in 1:length(all.combs[[j]])) {
            suggest[[s]] <- best.list[[all.combs[[j]][s]]][suggest.sample,]
          }
          
          suggest.c <- do.call(cbind, suggest)
        } else {
          suggest.c <- NULL
        }
        
        
        # Update GA.input

        # *-* Island GA ---------------------------------------------------------------
        
        if(isTRUE(GA.input_orig$gaisl)) {
          GA.inputs <- GA.prep(raster = raster.comb,
                               Results.dir = 'all.comb',
                               min.cat = GA.input_orig$inputs$min.cat,
                               max.cat = GA.input_orig$inputs$max.cat,
                               max.cont = GA.input_orig$inputs$max.cont,
                               cont.shape = NULL,
                               select.trans = s.trans,
                               method = GA.input_orig$inputs$method,
                               k.value = GA.input_orig$inputs$k.value,
                               pop.mult = GA.input_orig$inputs$pop.mult,
                               percent.elite = GA.input_orig$inputs$percent.elite,
                               type = GA.input_orig$inputs$type,
                               pcrossover = GA.input_orig$inputs$pcrossover,
                               pmutation = GA.input_orig$inputs$pmutation,
                               maxiter = GA.input_orig$inputs$maxiter,
                               run = GA.input_orig$inputs$run,
                               keepBest = GA.input_orig$inputs$keepBest,
                               population = GA.input_orig$inputs$population,
                               selection = GA.input_orig$inputs$selection,
                               crossover = GA.input_orig$inputs$crossover,
                               mutation = GA.input_orig$inputs$mutation,
                               pop.size = GA.input_orig$inputs$pop.size,
                               parallel = GA.input_orig$inputs$parallel,
                               gaisl = GA.input_orig$inputs$gaisl,
                               island.pop = GA.input_orig$inputs$island.pop,
                               numIslands = GA.input_orig$inputs$numIslands,
                               migrationRate = GA.input_orig$inputs$migrationRate,
                               migrationInterval = GA.input_orig$inputs$migrationInterval,
                               optim = GA.input_orig$inputs$optim,
                               optim.method = GA.input_orig$inputs$optim.method, 
                               poptim = GA.input_orig$inputs$poptim,
                               pressel = GA.input_orig$inputs$pressel,
                               control = GA.input_orig$inputs$control,
                               hessian = GA.input_orig$inputs$hessian,
                               seed = GA.input_orig$inputs$seed,
                               monitor = GA.input_orig$inputs$monitor,
                               quiet = GA.input_orig$inputs$quiet
          ) 
        } else { # *-* Standard GA ----------------------
          GA.inputs <- GA.prep(raster = raster.comb,
                               Results.dir = 'all.comb',
                               min.cat = GA.input_orig$inputs$min.cat,
                               max.cat = GA.input_orig$inputs$max.cat,
                               max.cont = GA.input_orig$inputs$max.cont,
                               cont.shape = NULL,
                               select.trans = s.trans,
                               method = GA.input_orig$inputs$method,
                               k.value = GA.input_orig$inputs$k.value,
                               pop.mult = GA.input_orig$inputs$pop.mult,
                               percent.elite = GA.input_orig$inputs$percent.elite,
                               type = GA.input_orig$inputs$type,
                               pcrossover = GA.input_orig$inputs$pcrossover,
                               pmutation = GA.input_orig$inputs$pmutation,
                               maxiter = GA.input_orig$inputs$maxiter,
                               run = GA.input_orig$inputs$run,
                               keepBest = GA.input_orig$inputs$keepBest,
                               population = GA.input_orig$inputs$population,
                               selection = GA.input_orig$inputs$selection,
                               crossover = GA.input_orig$inputs$crossover,
                               mutation = GA.input_orig$inputs$mutation,
                               pop.size = GA.input_orig$inputs$pop.size,
                               parallel = GA.input_orig$inputs$parallel,
                               optim = GA.input_orig$inputs$optim,
                               optim.method = GA.input_orig$inputs$optim.method, 
                               poptim = GA.input_orig$inputs$poptim,
                               pressel = GA.input_orig$inputs$pressel,
                               control = GA.input_orig$inputs$control,
                               hessian = GA.input_orig$inputs$hessian,
                               seed = GA.input_orig$inputs$seed,
                               monitor = GA.input_orig$inputs$monitor,
                               quiet = GA.input_orig$inputs$quiet
          ) 
        }
        
        ## Update suggests
        GA.inputs$SUGGESTS <- suggest.c
        
        # Update GA.input directories
        GA.inputs$Plots.dir <- paste0(results.dir,
                                      'rep_',i, 
                                      "/", 
                                      comb.names[[j]],
                                      "/")
        
        GA.inputs$Results.dir <- paste0(results.dir,
                                        'rep_',i, 
                                        "/", 
                                        comb.names[[j]],
                                        "/")
        
        if(!is.null(gdist.inputs)) {
          ms.results[[j]] <- MS_optim(gdist.inputs = gdist.inputs,
                                      GA.inputs = GA.inputs)
        } else if(!is.null(jl.inputs)) {
          ms.results[[j]] <- MS_optim(jl.inputs = jl.inputs,
                                      GA.inputs = GA.inputs)
        } else {
          stop("`gdist.inputs` or `jl.inputs` must be specified!!!")
        }
        
        all.cd[[(j + n_ss.cd)]] <- ms.results[[j]]$cd[[1]]
        ms.cd.names <- names(ms.results[[j]]$cd)
        ms.cd.name <- ""
        if (!is.null(ms.cd.names) && length(ms.cd.names) > 0) {
          ms.cd.name <- ms.cd.names[[1]]
        }
        if (is.null(ms.cd.name) || !nzchar(ms.cd.name) || ms.cd.name == "NULL") {
          ms.cd.name <- as.character(ms.results[[j]]$k$surface[[1]])
        }
        names(all.cd)[[j + n_ss.cd]] <- ms.cd.name
        
        ms.k[[j]] <- ms.results[[j]]$k
        
        AICc.tab_list[[j]] <- ms.results[[j]]$AICc.tab
    } # End combination loop
    
    # Convert combination lists to data frames
    if(is.null(ss.results)) {
      all.k <- plyr::ldply(ms.k)
      all.AICc <- plyr::ldply(AICc.tab_list)
    } else {
      all.k <- rbind(ss.results$k,
                     plyr::ldply(ms.k))
      all.AICc <- rbind(ss.results$AICc,
                        plyr::ldply(AICc.tab_list))
    }

    all.cd.names <- names(all.cd)
    if (is.null(all.cd.names)) {
      all.cd.names <- rep("", length(all.cd))
    }

    unnamed.cd <- which(is.na(all.cd.names) |
                          !nzchar(all.cd.names) |
                          all.cd.names == "NULL")
    if (length(unnamed.cd) > 0) {
      candidate.names <- as.character(all.k$surface)
      candidate.names <- candidate.names[nzchar(candidate.names) &
                                           candidate.names != "Null"]
      if (length(candidate.names) == length(all.cd)) {
        all.cd.names[unnamed.cd] <- candidate.names[unnamed.cd]
        names(all.cd) <- all.cd.names
      }
    }
    
    boot.k <- all.k[match(names(all.cd), all.k$surface), , drop = FALSE]

    if(anyNA(boot.k$surface) || nrow(boot.k) != length(all.cd)) {
      stop(
        "Failed to align optimized distance matrices with model metadata in ",
        "`all_comb()`. Please check single- and multi-surface outputs."
      )
    }
    
    all.AICc <- all.AICc %>% 
      dplyr::arrange(., AICc) %>%
      dplyr::mutate(., delta.AICc = AICc - min(AICc)) %>%
      dplyr::mutate(., weight = (exp(-0.5 * delta.AICc)) / sum(exp(-0.5 * delta.AICc))) %>%
      as.data.frame()
    
    
    # Bootstrap ----------------------------------------------------------
    if(!is.null(gdist.inputs)) {
      obs <- gdist.inputs$n.Pops
      genetic.mat <- matrix(0, obs, obs)
      genetic.mat[lower.tri(genetic.mat)] <- gdist.inputs$response
      
      boot.results <- Resist.boot(mod.names = boot.k[,1],
                                  dist.mat = all.cd,
                                  n.parameters = boot.k[,2],
                                  sample.prop = sample.prop,
                                  iters = iters,
                                  obs = obs,
                                  genetic.mat = genetic.mat)
    } else {
      obs <- jl.inputs$n.Pops
      genetic.mat <- matrix(0, obs, obs)
      genetic.mat[lower.tri(genetic.mat)] <- jl.inputs$response.all
      
      boot.results <- Resist.boot(mod.names = boot.k[,1],
                                  dist.mat = all.cd,
                                  n.parameters = boot.k[,2],
                                  sample.prop = sample.prop,
                                  iters = iters,
                                  obs = obs,
                                  genetic.mat = genetic.mat,
                                  keep = jl.inputs$keep)
    }
    
    # Write AICc table and Boot Results to replicate results directory
    write.table(x = all.AICc,
                paste0(results.dir,'rep_',i,"/",
                       "All_Combinations_Summary.csv"),
                row.names = F,
                col.names = T,
                sep = ",")
    
    write.table(x = as.data.frame(boot.results),
                paste0(results.dir,'rep_',i,"/",
                       "Bootstrap_Results.csv"),
                row.names = F,
                col.names = T,
                sep = ",")
    
    if(replicate > 1) {
      Results[[i]] <- list(summary.table = all.AICc,
                           boot.results = boot.results,
                           all.k = all.k,
                           all.cd = all.cd,
                           genetic.dist_mat = genetic.mat,
                           ss.results = ss.results,
                           ms.results = ms.results
      )
      names(Results)[i] <- paste0('rep_',i)
      
    } else {
      Results <- list(summary.table = all.AICc,
                      boot.results = boot.results,
                      all.k = all.k,
                      all.cd = all.cd,
                      genetic.dist_mat = genetic.mat,
                      ss.results = ss.results,
                      ms.results = ms.results
      )
    }
    
  } # Close replicate loop
  
  Results <- resga_add_class(Results, "resga_all_comb")
  return(Results)
  
} # End function
