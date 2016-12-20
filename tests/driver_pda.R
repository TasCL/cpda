## Use piecewise PDA -------------------------------------
library(pda); data(lba)
data.table::data.table(d)
str(plba)
# List of 4
# $ DT1 : num [1:695] 0.488 0.801 0.376 0.507 0.532 ...
# $ DT2 : num [1:305] 0.538 0.77 0.568 0.271 0.881 ...
# $ eDT1: num [1:7020] 0.475 0.346 0.42 0.401 0.368 ...
# $ eDT2: num [1:2980] 0.703 0.693 0.704 0.462 0.468 ...

## Use piecewise pda logLik to get likelihood for each data point
## This algorithm calculates gaussian kernel directly
## (1) First argument, plba$DT1, is the data.
## (2) Second argument, plba$eDT1, is the simulation.
## (3) The outputs are likelihoods (not log yet) for each data point
outputs <- logLik_pw(plba$DT1, plba$eDT1)

str(outputs)
sum(outputs)
sum(log(outputs))  ## 278.6095


## Use ordinary PDA -------------------------------------
## logLik_fft return summed, logged likelihood
logLik_fft(plba$DT1, plba$eDT1) ## 278.4902
logLik_fft(plba$DT2, plba$eDT2)

## Return all four outputs
## (1) summed, logged likelihood
## (2) likelihoods for each data point (ordinary PDA, not piecewise algorithm)
## (3) grid center (z)
## (4) PDA simulated histogram
tmp1 <- logLik_fft2(plba$DT2, plba$eDT2)
## manually piecewise using ordinary PDA
tmp2 <- logLik_fft2(plba$DT2[2], plba$eDT2)


for(i in 1:10) {
  tmp0 <- logLik_fft2(plba$DT2[i], plba$eDT2)
  cat("LL: ", round(tmp0$LL, 3),"\t")
  cat("PDF: ", round(tmp0$PDF, 3),"\n")
}


tmp2 <- density(plba$DT1)
tmp3 <- density(plba$DT1[1], bw=.02)
tmp3 <- density(plba$DT1[1])
tmp3$y

tmp0 <- logLik_fft(i, rnorm(1000))




## Plot normal
pVec <- c(mu=5, sigma=1)
y    <- sort(rnorm(1e5, 5, 1)) ## data
mcmcParams <- list(sigma=pVec[2], bandwidth=bwNRD0(y, 1), nmc=30,
                   ns=length(y), nchain=24, rp=.001, burnin=10, nthin=3,
                   start=1, gammaMult=2.38, ST=0, report=100)

setting <- c(mcmcParams$sigma,  mcmcParams$bandwidth, mcmcParams$ns,
             mcmcParams$nmc,    mcmcParams$nchain,    mcmcParams$rp,
             mcmcParams$burnin, mcmcParams$nthin,     mcmcParams$start,
             mcmcParams$gammaMult,   mcmcParams$ST, mcmcParams$report)

ll1 <- logLik_norm(y, pVec, setting)
ll2 <- logLik_norm2(y, pVec, setting)

str(ll2)

plot(ll2$z, ll2$PDF_hist, type="l", lty=2,
  main="Normal Distribution",xlab="x",ylab="L(x|beta)")
lines(y, ll2$PDF, col="red", lwd = 1.5)

##########################################################################


