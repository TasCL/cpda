rm(list=ls())
pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51, 
          muw2=3.69, t_delay=0.31, sv=1, swt=0.5)

data(lba)
dMat <- data.matrix(d); head(dMat)

set.seed(123)
samp   <- cpda::rplba(1e5, pVec); head(samp)
h   <- 0.8*bw.nrd0(samp[,2]); h
set.seed(123)
ll0 <- cpda::logLik_plba(dMat, pVec, 1e5, h=h, m=1); ll0

set.seed(123)
llList <- cpda::logLik_plba2(dMat, pVec, 1e5, h=h, m=1); 
str(llList)

## List of 4
## $ LL      : num -7277
## $ PDF     : num [1:1000, 1] -23 -23 -23 -23 -23 ...
## $ z       : num [1:2048, 1] 0.175 0.179 0.184 0.188 0.192 ...
## $ PDF_hist: num [1:2048, 1] 0 0 0 0 0 0 0 0 0 0 ...

## Here I demonstrate use rplba and logLik_fft to produce identical 
## result
set.seed(123)
DTMat  <- cpda::choiceDT(dMat, pVec)
time1  <- sort(DTMat[DTMat[,1] == 1, 2])
time2  <- sort(DTMat[DTMat[,1] == 2, 2])
time1_ <- sort(samp[samp[,1] == 1, 2])
time2_ <- sort(samp[samp[,1] == 2, 2])
ll1 <- cpda::logLik_fft(time1, time1_, h=h, m=1); ll1
ll2 <- cpda::logLik_fft(time2, time2_, h=h, m=1); ll2
print(ll1+ll2)
