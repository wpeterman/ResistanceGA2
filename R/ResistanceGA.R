#' @keywords internal
#' @details
#'
#' This package provides functions to optimize continuous and categorical
#' resistance surfaces using
#' \href{https://github.com/Circuitscape/Circuitscape.jl}{Circuitscape.jl}
#' (run through Julia) and Genetic Algorithms within R.
#' Resistance distances can alternatively be computed using the
#' \pkg{gdistance} package (least-cost paths / commute distances).
#'
#' Output includes: a summary table with AIC, AICc, conditional and marginal
#' R-squared values, and log-likelihood values for each optimized surface;
#' parameters that optimized each top model; coefficients from fitted mixed
#' effects models; response-curve plots; and diagnostic plots of model fit.
#'
#' @references
#' Peterman, W.E., G.M. Connette, R.D. Semlitsch, and L.S. Eggert. 2014.
#' Ecological resistance surfaces predict fine-scale genetic differentiation in
#' a terrestrial woodland salamander. Molecular Ecology 23:2402--2413.
#' \doi{10.1111/mec.12708}
#'
#' Peterman, W. E. 2018. ResistanceGA: An R package for the optimization of
#' resistance surfaces using genetic algorithms. Methods in Ecology and
#' Evolution 9:1638--1647. \doi{10.1111/2041-210X.12984}
#'
#' @author Bill Peterman \email{peterman.73@@osu.edu}
#'
#' @import GA ggplot2 gdistance terra
#' @importFrom lme4 mkMerMod lFormula glFormula mkGlmerDevfun optimizeGlmer mkLmerDevfun optimizeLmer GHrule
#' @importFrom utils combn menu
#' @importFrom ggExtra removeGrid ggMarginal
#' @importFrom Matrix fac2sparse drop0
#' @importFrom plyr arrange rbind.fill ldply create_progress_bar progress_text
#' @importFrom dplyr mutate group_by summarise filter tally left_join dense_rank
#' @importFrom MuMIn r.squaredGLMM
#' @importFrom spdep dnearneigh nb2mat
#' @importFrom grDevices dev.off tiff topo.colors
#' @importFrom graphics abline par
#' @importFrom stats AIC lm logLik qqline qqnorm resid residuals runif as.formula sigma fitted coef dist
#' @importFrom utils file_test read.csv read.delim read.table write.table
#' @importFrom JuliaConnectoR juliaEval juliaCall juliaImport
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %do% %dopar% foreach
"_PACKAGE"

# Quiet R CMD check notes for NSE variables
utils::globalVariables(c(".", "delta", "LL", "R2m", "RMSE", "desc",
                         "surface", "weight", "avg.rank", "n",
                         "cor.df", "pop1", "pop2", "AICc",
                         "delta.AICc"))


#' Data file. Simulated true resistance surfaces
#'
#' A \code{SpatRaster} with three layers generated from the \code{RandomFields}
#' package. These surfaces were used with individual-based genetic simulations
#' conducted using \code{PopGenReport}.
#'
#' \itemize{
#'   \item cat_true. 4-class categorical resistance surface. Reclass: 1=1, 2=75, 3=200, 4=25, 5=1.
#'   \item cont_true. Continuous resistance surface transformed via inverse monomolecular
#'     (shape = 2.75, max = 100).
#'   \item multi_true. Composite surface: sum of transformed categorical and continuous
#'     surfaces, rescaled to minimum of 1.
#' }
#'
#' @docType data
#' @name raster_true
#' @format \code{SpatRaster} with 3 layers
#' @usage data(raster_true)
#' @keywords datasets
#' 
NULL


#' Data file. Simulated original resistance surfaces
#'
#' A \code{SpatRaster} with two layers generated from the \code{RandomFields}
#' package. These are the raw (untransformed) surfaces used as optimization inputs.
#'
#' \itemize{
#'   \item cat_orig. 5-class categorical resistance surface.
#'   \item cont_orig. Continuous resistance surface.
#' }
#'
#' @docType data
#' @name raster_orig
#' @format \code{SpatRaster} with 2 layers
#' @usage data(raster_orig)
#' @keywords datasets
#' 
NULL


#' Data file. A list of three coordinate matrices
#'
#' 50 sample locations per list element corresponding to populations used in
#' genetic simulations. Each element is a two-column numeric matrix (x, y).
#'
#' \itemize{
#'   \item sample_cat. Sample locations for categorical surface simulation.
#'   \item sample_cont. Sample locations for continuous surface simulation.
#'   \item sample_multi. Sample locations for multi-surface simulation.
#' }
#'
#' @docType data
#' @name sample_pops
#' @format A list of length 3; each element is a 50 x 2 coordinate matrix.
#' @usage data(sample_pops)
#' @keywords datasets
NULL


#' Data file. A list of three matrices 
#' 
#' Each matrix depicts the pairwise genetic distance (measured as chord distance) between sample locations
#' 
#'  \itemize{    
#'    \item Dc_cat. Pairwise chord distance matrix simulated across the categorical resistance surface
#'    \item Dc_cont. Pairwise chord distance matrix simulated across the continuous resistance surface
#'    \item Dc_multi. Pairwise chord distance matrix simulated across the multivariate resistance surface
#'    }
#' 
#' @docType data
#' @name Dc_list
#' @format A list of length 3
#' @usage data(Dc_list)
#' @description  Sample file to be used with examples in the vignette
#' @keywords datasets
NULL


#' Data file. A list of three matrices 
#' 
#' Each matrix depicts the pairwise effective resistance distance (measured using `gdistance`) between sample locations across `resist_true` surfaces
#' 
#'  \itemize{      
#'    \item resist_cat. Pairwise effective distance matrix calculated across the categorical resistance surface
#'    \item resist_cont. Pairwise effective distance matrix calculated across the continuous resistance surface
#'    \item resist_cont. Pairwise effective distance matrix calculated across the multivariate resistance surface
#'    }
#' 
#' @docType data
#' @name resist_list
#' @format A list of length 3
#' @usage data(resist_list)
#' @description  Sample file to be used with examples in the vignette
#' @keywords datasets
NULL



#' Data file. Example sample location file
#' 
#' A data frame that can be saved as a .txt file for running examples in the vignette 
#' 
#' @docType data
#' @name samples
#' @format A 25 x 3 data frame
#' @usage data(samples)
#' @description  Sample file to be used with examples in the vignette
#' @keywords datasets
NULL

#' Data file. Example raster output from CIRCUITSCAPE Julia
#'
#' A \code{SpatRaster} produced by a Circuitscape Julia run.
#'
#' @docType data
#' @name jl_out
#' @format \code{SpatRaster} with 1 layer
#' @usage data(jl_out)
#' @keywords datasets
NULL

#' Data file. Simulated resistance surfaces
#'
#' A \code{SpatRaster} with three layers for use in vignette examples.
#'
#' \itemize{
#'   \item categorical. 3-class categorical resistance surface.
#'   \item continuous. Continuous resistance surface.
#'   \item feature. 2-class categorical (feature) resistance surface.
#' }
#'
#' @docType data
#' @name resistance_surfaces
#' @format \code{SpatRaster} with 3 layers
#' @usage data(resistance_surfaces)
#' @keywords datasets
#' 
NULL