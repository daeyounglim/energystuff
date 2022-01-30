## energystuff
This repository contains the R functions running the proposed algorithm in *A Hybrid Monitoring Procedure for Detecting Abnormality with Application to Energy Consumption Data (submitted)*. Functions for generating simulation data sets and reproducing the simulation results of the submitted paper are under the directory `R`; however, these simulation-related functions are not exported from the R package, nor are they documented.

Documentation for all other functions is available through the `help()` function—e.g., `help(hybridmonitor.multiple)`.

### Installation
+ Download the tar file `energystuff_1.0.tar.gz`.
+ Open R, set your working directory to where the tar file is located.
+ Run the following code to install `energystuff`:
```r
install.packages("energystuff_1.0.tar.gz", type="source", repo=NULL)
```

### Simulation-related functions
There are two files for reproducing the simulation results: `gendata_multiple.R` and `wrapper_multiple.R`. The suffix `_multiple` indicates the multiple-account version (see Lim et al. once published).

+ `gendata_multiple.R` - Program that generates an R `list` of `data.frame`s. Run `source(gendata_multiple.R)` to generate 10,000 simulation data sets each for size analysis and power analysis. The saved filenames will contain either `_size_` or `_power_` to indicate which analysis it is created for.
+ `wrapper_multiple.R` - This file contains an R function `simfun` that reads a suitable `RData` and runs the corresponding analysis—e.g., size analysis or power analysis. The function arguments are
    + `start` and `end` - the IDs of the first and last data sets for analysis. These arguments are intended to provide flexibility in choosing which data sets to include in the analysis. It also provides the ability to divvy up the data sets for distributed processing.
    + `seednumber=18007` - the seed number for random number generation used in the analysis.
    + `type="size"` - a character-type variable indicating the type of analysis, either `size` or `power`.
    + `save.file=TRUE` - a logical variable indicating whether to save the results as a file. If set to `FALSE`, the results will be remain as objects in the R environment.
    + `ncores=1` - the number of CPU cores used for parallel processing, provided that R is installed with OpenMP support.
    + `verbose=FALSE` - a logical variable indicating whether to print out the progress bar.

**Caution** - These simulations can take up to several hours even with a moderate number of data sets. Please take that in mind before running them to reproduce the results. These simulations have been distributed across 48 server nodes for our paper.