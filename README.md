# Probability Density Approximation

This implements probability density approximation, using Armadillo C++ library.

## Getting Started

```
require(cpda)

## Point-wise Example --------------------
library(cpda); data(lba)
## data.table::data.table(d)
str(plba)
## List of 4
## $ DT1 : num [1:695] 0.488 0.801 0.376 0.507 0.532 ...
## $ DT2 : num [1:305] 0.538 0.77 0.568 0.271 0.881 ...
## $ eDT1: num [1:7020] 0.475 0.346 0.42 0.401 0.368 ...
## $ eDT2: num [1:2980] 0.703 0.693 0.704 0.462 0.468 ...

## Use piecewise pda logLik to get likelihood for each data point
## This algorithm calculates gaussian kernel directly
## (1) First argument, plba$DT1, is the data.
## (2) Second argument, plba$eDT1, is the simulation.
## (3) The outputs are likelihoods (not log yet) for each data point

outputs <- logLik_pw(plba$DT1, plba$eDT1)
str(outputs)
sum(outputs)
sum(log(outputs))

## Gaussian Example --------------------
pVec <- c(mu=5, sigma=1)
y    <- sort(rnorm(1e5, pVec[1], pVec[2]))

ll1 <- logLik_norm(y, pVec, 1e5)
ll2 <- logLik_norm2(y, pVec, 1e5)
str(ll1); str(ll2)

plot(ll2$z, ll2$PDF_hist, type="l", lty=2,
main="Normal Distribution",xlab="x",ylab="L(x|beta)")
lines(y, ll2$PDF, col="red", lwd = 1.5)


## Piece-wise LBA Example -----------------------
rm(list=ls())
data(lba)
dMat <- data.matrix(d)
head(dMat)

pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51, muw2=3.69,
          t_delay=0.31, sv=1, swt=0.5)
## setting <- c(bandwidth=.02, ns=1e5, nmc=30, nchain=24, rp=.001, burnin=10,
## nthin=3, start=1, gammaMult=2.38, report=100)
tmp0 <- cpda::logLik_plba(dMat, pVec, 1e5); tmp0


```

## Installation 

```
## From github
devtools::install_github("TasCL/pda")
## From source: 
install.packages("pda_0.0.0.8.tar.gz", repos = NULL, type="source")
```

## Prerequisites
 - R (>= 3.0.2)
 - Rtools
 - Rcpp package (>= 0.12.3)
 
## References

* Holmes, W. (2015). A practical guide to the Probability Density
Approximation (PDA) with improved implementation and error characterization.
Journal of Mathematical Psychology, 68-69, 13--24,
doi: http://dx.doi.org/10.1016/j.jmp.2015.08.006.

