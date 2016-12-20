rm(list=ls())
cpda::logLik_pw(0, rnorm(1e6)) ## error
dnorm(0, log=TRUE)

h <- 0.8*bw.nrd0(rnorm(1e6)); h          
cpda::logLik_pw(0, rnorm(1e6), h=h) ## error
dnorm(0, log=TRUE)

cpda::logLik_pw(0, rnorm(1e6), h=h, m=1) ## error
dnorm(0, log=TRUE)

cpda::logLik_pw(0, rnorm(3e6)) ## error
dnorm(0, log=TRUE)

data(lba)
## data.table::data.table(d)
str(plba)

outputs <- logLik_pw(plba$DT1, plba$eDT1)
str(outputs)
sum(outputs)

rm(list=ls())
x    <- seq(-3, 3, length.out=100) ## Data
samp <- rnorm(1e7)              ## Monte Carlo simulation 
h <- 0.8*stats::bw.nrd0(samp); h          ## Define bin
system.time(pw1  <- logLik_pw(x, samp))
system.time(pw2  <- logLik_pw(x, samp, h=h, m=0.8))

system.time(pw3  <- logLik_pw(x, samp, h=h, m=1))
system.time(pw4  <- logLik_pw(x, samp, h=h, m=1, parallel=TRUE))

plot(x, pw1,type="l", lty="dotted")
lines(x, pw2, col="darkgreen", lty="dashed")
lines(x, pw3, col="blue", lty="dotdash")
lines(x, pw4, col="red", lty="longdash")
lines(x, dnorm(x, log=TRUE), col="red")


## First, I demonstrate how to use point-wise algorithm
## Note logLik_pw does logarithm internally.
system.time(pw1  <- logLik_pw(x, samp, h=h, m=1))
##   user  system elapsed 
##  3.480   0.120   3.605 

## Second, I demonstrate how to use KDE-FFT. logLik_fft2 retrieves 
## logged-likelihood for each data point  
system.time(fft1 <- logLik_fft2(x, samp, h=h, m=1)[["PDF"]])

## Third, I demonstrate how to build logLik_fft and normal random number
## generator together inside C++. Here I use arma::randn() function, 
## instead of R's Rf_rnorm. Enter ?logLik_norm2 to see help page
system.time(fft2 <- logLik_norm2(x, pVec=c(0,1), 1e6, h=h, m=1)[["PDF"]])


plot(x, pw1,type="l", lty="dotted")
lines(x, fft1, col="darkgreen", lty="dashed")
lines(x, fft2, col="blue", lty="dotdash")
lines(x, dnorm(x, log=TRUE), col="red")

data(lba)
logLik_fft(plba$DT1, plba$eDT1)
logLik_fft(plba$DT2, plba$eDT2)
plbaEg1 <- logLik_fft2(plba$DT1, plba$eDT1)
plbaEg2 <- logLik_fft2(plba$DT2, plba$eDT2)

str(plbaEg1)

system.time(pw1  <- logLik_pw(x, samp, h=h, m=1))
system.time(pw2  <- logLik_pw(x, samp, h=h, m=1))


plot(x, pw1,type="l", lty="dotted")
lines(x, pw2, col="green")
lines(x, dnorm(x, log=TRUE), col="red")


system.time(fft1 <- logLik_fft2(x, samp, h=h, m=1)[["PDF"]])
system.time(fft2 <- logLik_norm2(x, pVec=c(0,1), 1e6, h=h, m=1)[["PDF"]])

plot(x, pw1,type="l", lty="dotted")
lines(x, fft1, col="darkgreen", lty="dashed")
lines(x, fft2, col="blue", lty="dotdash")
lines(x, dnorm(x, log=TRUE), col="red")



system.time(pw2  <- logLik_pw(x, samp))
system.time(fft3 <- logLik_fft2(x, samp)[["PDF"]])
system.time(fft4 <- logLik_norm2(x, pVec=c(0,1), 1e6)[["PDF"]])

plot(x, pw2,type="l", lty="dotted")
lines(x, fft3, col="darkgreen", lty="dashed")
lines(x, fft4, col="blue", lty="dotdash")
lines(x, dnorm(x, log=TRUE), col="red")
