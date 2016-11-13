# Probability Density Approximation

This implements probability density approximation, using Armadillo C++ library.

## Getting Started

```
require(pda)

## Use piecewise LBA data as an example
data(lba)
bandwidth <- .02
nsample   <- 1e4
m <- min(c(plba$DT1, plba$DT2)) - 3 * bandwidth
M <- max(c(plba$DT1, plba$DT2)) + 3 * bandwidth
logLik_fft(plba$DT1, plba$eDT1, m, M, bandwidth, nsample)
logLik_fft(plba$DT2, plba$eDT2, m, M, bandwidth, nsample)
tmp1 <- logLik_fft2(plba$DT1, plba$eDT1, m, M, bandwidth, nsample)
tmp2 <- logLik_fft2(plba$DT2, plba$eDT2, m, M, bandwidth, nsample)
str(tmp1)

## List of 4
## $ LL      : num 33.9
## $ PDF     : num [1:695, 1] 1.761 0.445 1.368 1.681 1.542 ...
## $ z       : num [1:1024, 1] 0.157 0.159 0.16 0.161 0.162 ...
## $ PDF_hist: num [1:1025, 1] 0 0 0 0 0 0 0 0 0 0 ...
str(tmp2)
## List of 4
## $ LL      : num -323
## $ PDF     : num [1:305, 1] 0.5526 0.2489 0.5588 0.0579 0.3016 ...
## $ z       : num [1:1024, 1] 0.157 0.159 0.16 0.161 0.162 ...
## $ PDF_hist: num [1:1025, 1] 0 0 0 0 0 0 0 0 0 0 ...
```

## Installation 

```
## From github
devtools::install_github("TasCL/pda")
## From source: 
install.packages("pda_0.0.0.5.tar.gz", repos = NULL, type="source")
```

## Prerequisities
 - R (>= 3.0.2)
 - Rtools
 - Rcpp package (>= 0.12.3)
 
## References

* Holmes, W. (2015). A practical guide to the Probability Density
Approximation (PDA) with improved implementation and error characterization.
Journal of Mathematical Psychology, 68-69, 13--24,
doi: http://dx.doi.org/10.1016/j.jmp.2015.08.006.

