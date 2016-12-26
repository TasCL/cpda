rm(list=ls())
library(cpda)
data(lba)
str(d)
tmp <- data.matrix(d)
head(x)
str(x)
str(tmp)
pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, 
          t_ND=0.08, muw1=1.51, 
          muw2=3.69, t_delay=0.31, sv=1, swt=0.5)


dMat <- data.matrix(d); head(dMat)

set.seed(123)

set.seed(123)
ll0 <- cpda::logLik_plba2(dMat, pVec, 1e5);

x <- cbind(rep(0:1, each=100), rep(seq(.5, 2, length.out=100), 2))
fft3 <- logLik_plba2(x, pVec=pVec, ns=1e5)

head(x)
head(dMat)
