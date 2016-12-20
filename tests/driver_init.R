pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51,  
          muw2=3.69, t_delay=0.31,  sv=1, swt=0.5)

## Setting sequence is critical, too! 
setting <- c(bandwidth=.02, ns=1e5, nmc=30, nchain=24, rp=.001, burnin=10,
             nthin=3, start=1, gammaMult=2.38, report=100)
tmp0 <- init(pVec, setting)
str(tmp0)
tmp0$theta[,,1:2]
## theta is npar x nchain x nmc. Only 1st slice is filled
