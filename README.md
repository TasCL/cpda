# Probability Density Approximation

This package implements probability density approximation, using Armadillo C++  
and Open MP libraries.

## Getting Started

The following examples are retrieved from cpda help pages

```
require(cpda)

###################
## Example 1     ##
###################
## Simulate the standard normal case
n    <- 1e6
samp <- rnorm(n)

## This is the bandwidth calcuation used on the simulated samples for smoothing
h <- bw.nrd0(samp)


## Calcualte Log likelihoods
x <- seq(-3,3,length.out=100)

## By default the following does point wise (i.e., for each element in x) 
## smoothing using bandwidth = h*m, where m=0.8 by defualt as recommeded 
## by Holmes
system.time(pw <- cpda::logLik_pw(x,samp))
##    user  system elapsed 
##   2.488   0.136   2.623 

## NB1: Output is a vector
## NB2: Dont use paralell=TRUE, still being debugged!

## Here the smoothing is done using an fft based on all x, which makes things
## much faster as x gets bigger than logLik_pw(). This is useful for the case 
## where many RTs share the same parameter.
system.time(fft1 <- cpda::logLik_fft2(x,samp)[["PDF"]])
##    user  system elapsed
##   1.219   0.000   1.218

## NB1: Ouptut is a list (see ?logLik_fft2, logLik_fft returns just a single 
##      number for sum log like and to be used in fitting.
## NB2: [["PDF"]] return a 2 x n matrix, with column 1== data and 
##      column 2= PDF

# Nice close match (dont expect things much better than 1e-4)
sort(exp(pw)-exp(fft2[,2]))
##   [1] -1.524787e-04 -9.654363e-05 -9.020370e-05 -8.697464e-05 -8.358955e-05
##   [6] -8.210528e-05 -8.154130e-05 -7.789309e-05 -7.663537e-05 -6.890160e-05
## ...
##  [91]  6.069234e-05  6.333527e-05  6.368462e-05  7.328864e-05  7.723575e-05
##  [96]  8.236230e-05  9.316117e-05  9.365130e-05  1.463610e-04  1.494125e-04

## But everything plots as visually identical
png(filename="test_pw.png", width=1024, height=768)
plot(x,exp(pw),type="l", lty="longdash")
lines(x,exp(fft2[,2]),col="green", lty="dotted")
dev.off()

###################
## Example 2     ##
###################
## Here, I demonstrated how to use build-in routine binding logLik_fft with 
## rnorm to get Gaussian density. Enter ?logLik_norm2 to see help page
## pVec stands for parameter vector. For a Gaussian model, pVec stores mean
## and standard deviation.
system.time(fft2 <- cpda::logLik_norm2(x, pVec=c(0,1), n=n, h=h, m=0.8)[["PDF"]])
##   user  system elapsed 
##  1.304   0.040   1.349  

## NB: Again, [["PDF"]] return a 2 x n matrix, with column 1== data and 
##     column 2= PDF

## You should get all 4 lines overlaping on top of one another
## the tails may be a bit off.
plot(x, exp(pw),type="l", lty="longdash")
lines(x, exp(fft1[,2]), col="darkgreen", lty="dashed")
lines(x, exp(fft2[,2]), col="blue", lty="dotdash")
lines(x, dnorm(x, log=FALSE), col="red")

###################
## Example 3     ##
###################
## Use piecewise LBA data as an example
data(lba)
logLik_fft(plba$DT1, plba$eDT1)
logLik_fft(plba$DT2, plba$eDT2)
plbaEg1 <- logLik_fft2(plba$DT1, plba$eDT1)
plbaEg2 <- logLik_fft2(plba$DT2, plba$eDT2)

str(plbaEg1)
List of 4
## $ LL      : num 280
## $ PDF     : num [1:695, 1:2] 0.217 0.223 0.228 0.23 0.233 ...
## $ z       : num [1:1024, 1] 0.128 0.129 0.13 0.131 0.133 ...
## $ PDF_hist: num [1:1024, 1] 0 0 0 0 0 0 0 0 0 0 ...
str(plbaEg2)
## List of 4
##  $ LL      : num 45.3
##  $ PDF     : num [1:305, 1:2] 0.271 0.311 0.317 0.318 0.328 ...
##  $ z       : num [1:1024, 1] 0.121 0.124 0.126 0.129 0.131 ...
##  $ PDF_hist: num [1:1024, 1] 0 0 0 0 0 0 0 0 0 0 ...

###################
## Example 4     ##
###################
## Gaussian distribution
pVec <- c(mu=5, sigma=1)
y    <- sort(rnorm(1e5, pVec[1], pVec[2]))
ll1  <- logLik_norm(y, pVec, 1e5)
ll2  <- logLik_norm2(y, pVec, 1e5)
str(ll1); str(ll2)
## num -142255
## List of 4
##  $ LL      : num -142241
##  $ PDF     : num [1:100000, 1:2] 0.74 0.838 0.879 1.151 1.162 ...
##  $ z       : num [1:1024, 1] 0.503 0.512 0.521 0.529 0.538 ...
##  $ PDF_hist: num [1:1024, 1] 0 0 0 0 0 0 0 0 0 0 ...
 
plot(ll2$z, ll2$PDF_hist, type="l", lty=2,
main="Normal Distribution",xlab="x",ylab="L(x|beta)")
lines(y, exp(ll2$PDF[,2]), col="red", lwd = 1.5)


###################
## Example 5     ##
###################
## LBA Example 
rm(list=ls())
## Demo identical outputs from manually building rlba and logLik_fft, 
## comparing to logLik_lba 

## First retrieve example data set
data(lba)  
d$R <- ifelse(d$Response==0, 1, 2)  ## convert 0 & 1 accumulator to 1 & 2
dMat <- data.matrix(data.frame(R=d$R, RT=d$ResponseTime))
head(dMat)
##      R        RT
## [1,] 1 0.5676343
## [2,] 2 0.6183177
## [3,] 1 0.8806298
## [4,] 1 0.4563023
## [5,] 2 0.8496136
## [6,] 1 0.5866219

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
install.packages("pda_0.0.1.5.tar.gz", repos = NULL, type="source")
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

