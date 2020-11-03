# tantale package
## An integrated collection of functions for TALE minning and analysis

<p align="center">
  <img src="./man/figures/pipeline.svg">
</p>



- A TALE-oriented OOP framework:
    - A TALE class and associated methods (to be done)


- TALE mining in bacterial sequences:
    - Wrapper around annotale
    - tellTale
    - Analysis tools for RVD inventory, repeat lenght


- TALEs classification, phylogeny:
    - Wrappers around distal, functal, annotale
    - TALE groups inference
    - Multiple alignments plotting  


- TALE targets mining:
    - Wrappers around target predictors
    - General parser for results aggregation
    - Connector with daTALbase

**NOTE** :  
    - create vignette by `devtools:install(build_vignettes = TRUE)`  
    - the perl lib is not included in the installed package, so `runDistal` doesn't work ... but I don't understand why `talvez()` works when you didn't specify a perl lib??

    (https://github.com/marcschwartz/WriteXLS,
      https://github.com/cran/gdata/blob/master/R/installXLSXsupport.R)
