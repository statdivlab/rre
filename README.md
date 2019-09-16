# rre 
Replicate Richness Estimation

**Authors:** Alex Paynter and [Amy Willis](http://statisticaldiversitylab.com/team/amy-willis)

## Introduction

`rre` is an R package for estimating species richness where multiple frequency count tables have been collected.  All methods in this package assume a negative binomial model and explore the use of penalization to leverage the replicate data for improved estimation. [`rre_sims`](https://github.com/statdivlab/rre_sims) is a related package containing code for simulations and data analyses.

## Installation

`rre` can be installed from GitHub by running the following code in an R console:

```r
# install.packages("devtools") # run if devtools is not already installed.
devtools::install_github("statdivlab/rre")
library(rre)
```

## Example

In this example we will demonstrate the use of two `rre` methods using simulated data.  The two methods we show are the negative binomial MLE without regularization and a method which uses goodness of fit to tune a penalization term.

### Generate data
```r
library(rre)

# Generate 2 replicates from a population with 100 species, 
#   with the negative binomial parameters alpha = 0.1 and delta = 0.1.
set.seed(12)
list_of_fct <- nb_fct_simulation(ccc = 100, alpha = 0.1, delta = 0.1, r = 2)
```

### MLE
```r
# The function unregularized_mle optimizes the negative binomial 
# likelihood with no penalization:
mle_result <- unregularized_mle(list_of_fct)
# To extract the MLE we extract the 'best' element:
mle_result$best # ccc is the richness estimate, in this case 79
```

### Goodness of fit method
```r
# gof_criterion is one of our novel methods, which uses penalization tuned 
#   by a goodness of fit metric to estimate richness:
gof_result <- gof_criterion(list_of_fct)
# Extract the 'best' element:
gof_result$best # ccc is the richness estimate, in this case 79
```

## Issues

Our penalization tuning methods use grids for both richness and the penalization parameter.  A warning message is printed when the estimate is near the upper end of either of these grids, with a recommended course of action.

If you encounter issues please let us know by [filing an issue](https:://github.com/statdivlab/rre/issues).