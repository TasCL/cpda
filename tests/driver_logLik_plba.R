rm(list=ls())
data(lba)
dMat <- data.matrix(d)
head(dMat)

pVec <- c(A=1.51, muv1=3.32, muw1=1.51, muv2=2.24, muw2=3.69,
          t_delay=0.31, t_ND=0.08, b=2.7, sv=1, swt=0.5)
tmp0 <- cpda::logLik_plba(dMat, pVec, 1e6); tmp0


rm(list=ls())
data(lba)
dMat <- data.matrix(d)
head(dMat)

pVec <- c(b1=1, b2=1, A1=.5, A2=.5, mu1=2.4, mu2=1.6, sigma1=1, sigma2=1.2,
          t01=.5, t02=.5)
tmp0 <- cpda::logLik_lba(dMat, pVec, 1e8); tmp0


## Testing codes
vec <- 1:10
logLik_fft(1, rnorm(1e6))
logLik_fft(1, rnorm(1e6), h=0.08)
logLik_fft(1, rnorm(1e6), h=0.08, m=0.7)

tmp <- logLik_fft2(vec, rnorm(2e7)); log(tmp$PDF)
log(dnorm(vec))

for(i in 1:length(vec)) {
  tmp <- logLik_fft(vec[i], rnorm(1e7), h=0.008); 
  cat(tmp,"\n")
}

tmp <- logLik_fft2(vec, rnorm(1e7), h=0.008); 
log(tmp$PDF)



x  <- rnorm(1e3)
samp <- rnorm(1e6)
h <- 0.8*bw.nrd0(samp)
log(logLik_pw(1, samp, h))
logLik_fft(1, samp, h)


## add bandwidth and multiplier 

fft <- exp(logLik_fft2(x, samp)[["PDF"]])
pw <- exp(logLik_pw())
plot(x, fft)

log(dnorm(1))


pVec <- c(mu=5, sigma=1)
y    <- sort(rnorm(1e3, pVec[1], pVec[2]))

ll1 <- logLik_norm(y, pVec, 1e6)
ll2 <- logLik_norm2(y, pVec, 1e6)
str(ll1); str(ll2)

plot(ll2$z, ll2$PDF_hist, type="l", lty=2,
     main="Normal Distribution",xlab="x",ylab="L(x|beta)")
lines(y, ll2$PDF, col="red", lwd = 1.5)

