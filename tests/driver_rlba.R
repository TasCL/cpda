## rlba --------
n    <- 10
pVec <- c(b1=1, b2=1, A1=.5, A2=.5, mu1=2.4, mu2=1.6, sigma1=1, sigma2=1.2,
          t01=.5, t02=.5)

dat1 <- cpda::rlba(n, pVec)
dat2 <- rtdists::rLBA(n, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2))
head(dat1)
head(dat2)

h <-cpda::bwNRD0(plba$DT1, 0.8)

## rplba --------
rm(list=ls())
data(lba)
dMat <- data.matrix(d)
head(dMat)

pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51, muw2=3.69,
          t_delay=0.31, sv=1, swt=0.5)

tmp0 <- cpda::logLik_plba(dMat, pVec, 1e5); tmp0


pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51, 
          muw2=3.69, t_delay=0.31, sv=1, swt=0.5)
DTs <- cpda::rplba(n=1e4, pVec)
head(DTs)



## choiceDT --------
rm(list=ls())
data(lba)
pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51,  
          muw2=3.69, t_delay=0.31,  sv=1, swt=0.5)

DTs <- cpda::choiceDT(data.matrix(d), pVec)
head(DTs)

## Wrong one
pVecWrong <- c(A=1.51, muv1=3.32, muw1=1.51, muv2=2.24, muw2=3.69,
          t_delay=0.31, t_ND=0.08, b=2.7, sv=1, swt=0.5)
DTX <- cpda::choiceDT(data.matrix(d), pVecWrong)
head(DTX)


