# ResistanceGA2

### An R package to optimize resistance surfaces using Genetic Algorithms. Both continuous and categorical surfaces can be optimized using these functions. Additionally, it is possible to simultaneously optimize multiple resistance surfaces at the same time to generate novel resistance surfaces. Resistance distances can be calculated as cost distances (least cost path) between points, or as circuit-based resistance distances calculated using CIRCUITSCAPE.

To install this package, execute the following commands in R:

```         
# Install 'devtools' package, if needed
if(!("devtools" %in% list.files(.libPaths()))) {
    install.packages("devtools", repo = "http://cran.rstudio.com", dep = TRUE) 
} 

devtools::install_github("wpeterman/ResistanceGA2", build_vignettes = TRUE) # Download package

library(ResistanceGA2) # Loads package and the other dependnecies
```

Once the package is installed, you can view the 'Getting Started' vignette in R.

------------------------------------------------------------------------

### Other notes

-   This is large and ongoing update. `terra` is now fully (and only) supported. There is also a switch away from `JuliaCall` to `JuliaConnectoR`. Please make me aware of issues / errors as you encounter them.

-   Public spatial inputs are now terra-only: resistance surfaces should be supplied as `terra::SpatRaster` objects and sample locations as `terra::SpatVector` point layers. `gdistance` support is retained, with internal conversion to the raster classes that package still requires.

-   The standalone Circuitscape `.exe` workflow is no longer supported. Circuitscape runs are handled through Julia and `JuliaConnectoR`.

-   For continuous surfaces, it is strongly recommended that transformation searches be restricted to Monomolecular families unless there is a strong biological or ecological reason to expect a unimodal or quadratic-type relationship. Exploring Ricker families more broadly can increase the risk of overfitting.

-   It is highly recommended that you install Julia and the CIRCUITSCAPE Julia package. General instructions [**here**](https://petermanresearch.weebly.com/uploads/2/5/9/2/25926970/julia_guide.pdf "Julia Guide"). There is also a `Julia_Guide` vignette with the package now.

------------------------------------------------------------------------

This approach has been developed from the methods first utilized in Peterman et al. (2014). The first formal analysis using ResistanceGA was Ruiz-López et al. (2016). The primary citation for the package is Peterman (2018), Methods in Ecology.

Peterman, W.E., G.M. Connette, R.D. Semlitsch, and L.S. Eggert. 2014. Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology 23:2402–2413. [**PDF**](http://petermanresearch.weebly.com/uploads/2/5/9/2/25926970/peterman_et_al._2014--mec.pdf "Peterman et al.")

Peterman, W. E. 2018. ResistanceGA: An R package for the optimization of resistance surfaces using genetic algorithms. Methods in Ecology and Evolution 9, 1638–1647. <doi:10.1111/2041-210X.12984>. [**PDF**](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.12984 "MEE Publication")

Ruiz-López, M.J., Barelli, C., Rovero, F., Hodges, K., Roos, C., Peterman, W.E., Ting, N., 2016. A novel landscape genetic approach demonstrates the effects of human disturbance on the Udzungwa red colobus monkey (*Procolobus gordonorum*). Heredity, Ruiz-López 116, 167–176. <https://doi.org/10.1038/hdy.2015.82>
