rm(list=ls())



## Here I demonstrate use rplba and logLik_fft to produce identical 
## result
set.seed(123)
samp <- cpda::rlba(1e5, pVec); head(samp)
h    <- 0.8*bw.nrd0(samp[,2]); h

set.seed(123)
tmp0 <- cpda::logLik_lba(dMat, pVec, 1e5, h, 1); tmp0


dt1  <- sort(dMat[dMat[,1] == 1, 2]) - pVec[9]
dt2  <- sort(dMat[dMat[,1] == 2, 2]) - pVec[10]
dt1_ <- sort(samp[samp[,1] == 1, 2]) - pVec[9]
dt2_ <- sort(samp[samp[,1] == 2, 2]) - pVec[10]

ll1 <- cpda::logLik_fft(dt1, dt1_, h, 1); ll1
ll2 <- cpda::logLik_fft(dt2, dt2_, h, 1); ll2
print(ll1+ll2)



tmp1 <- cpda::logLik_lba2(dMat, pVec, 1e5); tmp1
str(tmp1)


dim(dMat)

