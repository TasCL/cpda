# Probability Density Approximation

This is the sister package of _ppda_ for conducting probability density 
approximation.  cpda provides separate modules for fast random number 
generation and likelihood functions, for example ex-Gaussian distribution. This
design allows the user who are not familar with C/C++ to build their PDA model 
in R. To gain the benefit of massive MCMC, we recommend using GPU computation.  


## Installation 

```
## From github
devtools::install_github("TasCL/cpda")
## From source: 
install.packages("pda_0.1.1.tar.gz", repos = NULL, type="source")

```

## Prerequisites
 - R (>= 3.0.2)
 - Rtools
 - Rcpp package (>= 0.12.3)
 - RcppArmadillo (>= 0.6.700.6.0)
 
## Contributors

- Yi-Shin Lin <yishin.lin@utas.edu.au> 
- Andrew Heathcote 
- William Holmes 

## References

* Holmes, W. (2015). A practical guide to the Probability Density
Approximation (PDA) with improved implementation and error characterization.
Journal of Mathematical Psychology, 68-69, 13--24,
doi: http://dx.doi.org/10.1016/j.jmp.2015.08.006.

