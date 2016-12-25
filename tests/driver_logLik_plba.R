###################
## Example 1     ##
###################
## I demonstrate how to use the build-in logLik_plba to calculate
## pLBA densities. set.seed(123) is to produce the same result.
rm(list=ls())
pVec <- c(A1=1.51, A2=1.51, b1=2.7, b2=2.7, v1=3.32, v2=2.24,
          w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.31, 
          swt=0.5, t0=0.08)

data(lba)
d$Response <- ifelse(d$Response==0, 1, 2)
dMat <- data.matrix(d); head(dMat)
n <- 1e5

set.seed(123)
ll0 <- cpda::logLik_plba(dMat, pVec, 1e5, h=h, m=1); ll0

set.seed(123)
llList <- cpda::logLik_plba2(dMat, pVec, 1e5, h=h, m=1); str(llList)
# List of 4
# $ LL      : num 327
# $ PDF     : num [1:1000, 1:4] 1 1 1 1 1 1 1 1 1 1 ...
# $ z       : num [1:2048, 1] 0.13 0.135 0.14 0.146 0.151 ...
# $ PDF_hist: num [1:2048, 1] 0 0 0 0 0 ...

###################
## Example 2     ##
###################
## Secondly, I demonstrate use rplba2 and logLik_fft to produce identical 
## result
rm(list=ls())
x <- cbind(rep(1:2,each=100),rep(seq(.5,2,length.out=100),2))
pVec <- c(A1=1.51, A2=1.51, b1=2.7, b2=2.7, v1=3.32, v2=2.24,
          w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.31, 
          swt=0.5, t0=0.08)
dt1 <- x[x[,1]==1,2] - pVec[15]
dt2 <- x[x[,1]==2,2] - pVec[15]

set.seed(123)
n <- 1e5
samp <- cpda::rplba2(n, pVec); head(samp)
dt1_ <- samp[samp[,1]==1,2] - pVec[15]
dt2_ <- samp[samp[,1]==2,2] - pVec[15]
logLik_fft(dt1, dt1_) + logLik_fft(dt2, dt2_)

set.seed(123)
cpda::logLik_plba(x, pVec)

fft1 <- cpda::logLik_plba2(x, pVec)
str(fft1)
head(fft1$PDF)
## choice RT PDF
tmp0 <- fft1$PDF[fft1$PDF[,1] == 1, 2]
tmp1 <- sort(x[x[,1]==1,2])
all(tmp0==tmp1)

tmp3 <- fft1$PDF[fft1$PDF[,1] == 2,2]
tmp4 <- sort(x[x[,1]==1,2])
all(tmp3==tmp4)



