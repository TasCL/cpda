pVec <- c(mu=5, sigma=1)
y    <- sort(rnorm(1e5, pVec[1], pVec[2]))

ll1 <- logLik_norm(y, pVec, 1e5)
ll2 <- logLik_norm2(y, pVec, 1e5)
str(ll1); str(ll2)

plot(ll2$z, ll2$PDF_hist, type="l", lty=2,
     main="Normal Distribution",xlab="x",ylab="L(x|beta)")
lines(y, exp(ll2$PDF), col="red", lwd = 1.5)
