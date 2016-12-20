library(cpda)
pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51,  
           muw2=3.69, t_delay=0.31,  sv=1, swt=0.5)
setting <- c(bandwidth=.02, ns=1e4, nmc=6, nchain=24, rp=.001, burnin=10,
             nthin=3, start=1, gammaMult=2.38, report=20)


samples <- init(pVec, setting)
useTheta <- samples$theta[,,1]

str(samples)
samples$opVec
data(lba)
dMat <- data.matrix(d)

tmp0 <- run(dMat, pVec, samples, setting)
tmp0
