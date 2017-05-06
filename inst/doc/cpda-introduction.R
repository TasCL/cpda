### R code from vignette source 'cpda-introduction.Rnw'

###################################################
### code chunk number 1: cpda-introduction.Rnw:
###################################################
library(cpda)
set.seed(123)
n   <- 1e5
obs <- rnorm(n)
h   <- bw.nrd0(obs);              ## KDE bandwidth
x   <- seq(-3, 3, length.out=100) ## Support

system.time(fft1 <- logLik_pw(x, obs))
##  user  system elapsed
## 0.328   0.000   0.323


###################################################
### code chunk number 2: cpda-introduction.Rnw:
###################################################
system.time(fft2 <- logLik_fft2(x, obs)[["PDF"]])
##  user  system elapsed
## 0.124   0.000   0.122

head(fft2)
##           [,1]      [,2]
## [1,] -3.000000 -5.420678
## [2,] -2.939394 -5.240477
## [3,] -2.878788 -5.048910
## [4,] -2.818182 -4.850126
## [5,] -2.757576 -4.671268
## [6,] -2.696970 -4.498645


###################################################
### code chunk number 3: cpda-introduction.Rnw:
###################################################
png("gaussian.png")
plot(x,  exp(fft1), type="l", lty="dotted", xlab="Observation", ylab="Density")
lines(x, exp(fft2[,2]), col="green", lty="dotdash")
lines(x, dnorm(x),      col="red",   lty="longdash")
dev.off()
