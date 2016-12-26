# Probability Density Approximation

This package implements probability density approximation, using Armadillo C++ and 
Open MP libraries.

## Getting Started

The following examples are retrieved from cpda help pages

```
require(cpda)

###################
## Example 1     ##
###################
x <- seq(-3, 3, length.out=100) ## Data
samp <- rnorm(1e6)              ## Monte Carlo simulation 
h <- 0.8*bw.nrd0(samp)          ## Define bandwidth using R's bw.nrd0

## First, I demonstrate how to use point-wise logLik
## Note logLik_pw returns log-likelihood.
system.time(pw1  <- logLik_pw(x, samp, h=h, m=1))
##   user  system elapsed 
##  3.480   0.120   3.605 

## Second, I demonstrate how to use KDE-FFT. logLik_fft2 retrieves 
## log-likelihoods for all data point  
system.time(fft1 <- logLik_fft2(x, samp, h=h, m=1)[["PDF"]])
##   user  system elapsed 
##  1.056   0.024   1.079 

## Third, I demonstrate how to build-in routine binding logLik_fft and 
## rnorm to get Gaussian density. Enter ?logLik_norm2 to see help page
## pVec stands for parameter vector. For a Gaussian model, pVec stores mean
## and standard deviation.
system.time(fft2 <- logLik_norm2(x, pVec=c(0,1), 1e6, h=h, m=1)[["PDF"]])
##   user  system elapsed 
##  1.304   0.040   1.349  

## You should get all 4 lines overlaping one another
## the tails may be a bit off.
plot(x, pw1,type="l", lty="dotted")
lines(x, fft1, col="darkgreen", lty="dashed")
lines(x, fft2, col="blue", lty="dotdash")
lines(x, dnorm(x, log=TRUE), col="red")

###################
## Example 2     ##
###################
## Use piecewise LBA data as an example
data(lba)
logLik_fft(plba$DT1, plba$eDT1)
logLik_fft(plba$DT2, plba$eDT2)
plbaEg1 <- logLik_fft2(plba$DT1, plba$eDT1)
plbaEg2 <- logLik_fft2(plba$DT2, plba$eDT2)

str(plbaEg1)
## List of 4
## $ LL      : num 280
## $ PDF     : num [1:695, 1] 0.918 -0.456 0.668 0.872 0.788 ...
## $ z       : num [1:1024, 1] 0.128 0.129 0.13 0.131 0.133 ...
## $ PDF_hist: num [1:1024, 1] 0 0 0 0 0 0 0 0 0 0 ...
str(plbaEg2)
## List of 4
## $ LL      : num 45.3
## $ PDF     : num [1:305, 1] 0.6396 -0.1174 0.6178 -1.3723 -0.0273 ...
## $ z       : num [1:1024, 1] 0.121 0.124 0.126 0.129 0.131 ...
## $ PDF_hist: num [1:1024, 1] 0 0 0 0 0 0 0 0 0 0 ...

## Gaussian Example --------------------
pVec <- c(mu=5, sigma=1)
y    <- sort(rnorm(1e5, pVec[1], pVec[2]))

ll1 <- logLik_norm(y, pVec, 1e5)
ll2 <- logLik_norm2(y, pVec, 1e5)
str(ll1); str(ll2)

plot(ll2$z, ll2$PDF_hist, type="l", lty=2,
main="Normal Distribution",xlab="x",ylab="L(x|beta)")
lines(y, exp(ll2$PDF), col="red", lwd = 1.5)


## LBA Example -----------------------
rm(list=ls())
## Demo identical outputs from manually building rlba and logLik_fft, 
## comparing to logLik_lba 

## First retrieve example data set
data(lba)  
d$R <- ifelse(d$Response==0, 1, 2)  ## convert 0 & 1 accumulator to 1 & 2
dMat <- data.matrix(data.frame(R=d$R, RT=d$ResponseTime))
head(dMat)
 
## LBA parameter vector. The sequence is critical. 
pVec <- c(b1=1, b2=1, A1=.5, A2=.5, mu1=2.4, mu2=1.6, sigma1=1, sigma2=1.2,
          t01=.5, t02=.5)
           
set.seed(123)  ## make sure using identical simulations
samp <- cpda::rlba(1e5, pVec); head(samp)
h    <- 0.8*bw.nrd0(samp[,2]); h

## logLik_lba simulates internally, so set.seed ahead to make sure using
## identical simulations
set.seed(123)  
tmp0 <- cpda::logLik_lba(dMat, pVec, 1e5, h, 1); tmp0 ## -3496.88

## Manually calculate empirical DTs for accumulator 1 and accumualtor 2
dt1  <- sort(dMat[dMat[,1] == 1, 2]) - pVec[9]
dt2  <- sort(dMat[dMat[,1] == 2, 2]) - pVec[10]
dt1_ <- sort(samp[samp[,1] == 1, 2]) - pVec[9]
dt2_ <- sort(samp[samp[,1] == 2, 2]) - pVec[10]
    
## Calculating log-likelihoods separately for each accumulator
ll1 <- cpda::logLik_fft(dt1, dt1_, h, 1); ll1
ll2 <- cpda::logLik_fft(dt2, dt2_, h, 1); ll2
print(ll1+ll2) ## -3496.88

```

## Installation 

```
## From github
devtools::install_github("TasCL/pda")
## From source: 
install.packages("pda_0.0.1.4.tar.gz", repos = NULL, type="source")
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

