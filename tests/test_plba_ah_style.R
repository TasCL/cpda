## Test codes for logLik_pw and logLik_fft2 ------------------------------------
## Simulate the standard normal case
n    <- 1e6
samp <- rnorm(n)

## This is the bandwidth calcuation used on the simulated samples for smoothing
h <- bw.nrd0(samp)

## Calcualte Log likelihoods
x <- seq(-3,3, length.out=100)

## By default the following does point wise (i.e., for each x) smoothing using
## bandwidth = h*m, where m=0.8 by defualt as recommeded by Holmes
system.time(pw <- cpda::logLik_pw(x,samp))
##    user  system elapsed
##   3.672   0.000   3.676
##   2.488   0.136   2.623 
## NB1: Output is a vector
## NB2: Do not use paralell=TRUE, still being debugged!
str(pw)
## num [1:100] -5.47 -5.27 -5.06 -4.89 -4.71 ...

## See R help page. 
?cpda::logLik_pw
## The user can change h and m too, if s/he does not wish to use a default 
## bandwidth. 
system.time(pw_tmp <- cpda::logLik_pw(x, samp, h=h, m=0.6))
## 2.160   0.136   2.295 
str(pw_tmp)
## num [1:100] -5.41 -5.22 -5.02 -4.86 -4.73 ...


## Here the smoothing is done using an fft-based method on all x, which makes 
## things much faster as x gets bigger. This is useful for the case where many 
## RTs share the same parameter.
system.time(fft2 <- cpda::logLik_fft2(x,samp)[["PDF"]])
##  user  system elapsed
## 1.219   0.000   1.218
## 0.760   0.012   0.773 

head(fft2)
##           [,1]      [,2]
## [1,] -3.000000 -5.403123
## [2,] -2.939394 -5.218121

## The user can request high grid precision too; in this case, it does not 
## make significant difference.
system.time(fft2.2 <- cpda::logLik_fft2(x,samp, p=12)[["PDF"]])
## 2.064   0.008   2.071 
head(fft2.2)
##           [,1]      [,2]
## [1,] -3.000000 -5.404047
## [2,] -2.939394 -5.215409

## See help page for logLik_fft2
?logLik_fft2
## NB1: Ouptut is a list (see ?logLik_fft2, logLik_fft returns just a single
## number for sum log like and to be used in fitting.)
## NB2: ["PDF"] return a 2 x n matrix, with column 1 is data and column 2 is PDF

# Nice close match (dont expect things much better than 1e-4)
sort(exp(pw)-exp(fft2[,2]))
##   [1] -1.524787e-04 -9.654363e-05 -9.020370e-05 -8.697464e-05 -8.358955e-05
##   [6] -8.210528e-05 -8.154130e-05 -7.789309e-05 -7.663537e-05 -6.890160e-05
## ...
##  [91]  6.069234e-05  6.333527e-05  6.368462e-05  7.328864e-05  7.723575e-05
##  [96]  8.236230e-05  9.316117e-05  9.365130e-05  1.463610e-04  1.494125e-04

## But everything plots as visually identical
png(filename="test_pw.png", width=1024, height=768)
plot(x,exp(pw),type="l", lty="longdash")
lines(x,exp(fft2[,2]),col="green", lty="dotted")
dev.off()

## Test codes for rplba---------------------------------------------------------

#################20
## Example 1     20
#################20
## Set up three parameter vectors testing 3 conditions:
## 1. rD (rate delay) == tD (threshold delay)
## 2. rD (rate delay) < tD (threshold delay)
## 3. rD (rate delay) > tD (threshold delay)
pVec3.1 <- c(A1=1.51, A2=1.51, B1=1.2, B2=1.2, C1=.3, C2=.3, v1=3.32, 
             v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.1,
             tD=.1, swt=0.5, t0=0.08)
pVec3.2 <- c(A1=1.51, A2=1.51, B1=1.2, B2=1.2, C1=.3, C2=.3, v1=3.32, 
             v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.1,
             tD=.15, swt=0.5, t0=0.08)
pVec3.3 <- c(A1=1.51, A2=1.51, B1=1.2, B2=1.2, C1=.3, C2=.3, v1=3.32, 
             v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.15,
             tD=.1, swt=0.5, t0=0.08)

n <- 1e5 ## This is enough to show overlapping densities
set.seed(123); system.time(dat5.1 <- cpda::rplbaR(n, pVec3.1))  ## R version
set.seed(123); system.time(dat5.2 <- cpda::rplbaR(n, pVec3.2))  
set.seed(123); system.time(dat5.3 <- cpda::rplbaR(n, pVec3.3))  
set.seed(123); system.time(dat6.1 <- cpda::rplba( n, pVec3.1))  ## C++ version
set.seed(123); system.time(dat6.2 <- cpda::rplba( n, pVec3.2))
set.seed(123); system.time(dat6.3 <- cpda::rplba( n, pVec3.3))

## Around 8-9 times difference in speed
## set.seed(123); system.time(dat5.1 <- cpda::rplbaR(n, pVec3.1))
##    user  system elapsed 
##   0.216   0.000   0.216 
## set.seed(123); system.time(dat5.2 <- cpda::rplbaR(n, pVec3.2))
##    user  system elapsed 
##   0.216   0.000   0.216 
## set.seed(123); system.time(dat5.3 <- cpda::rplbaR(n, pVec3.3))
##    user  system elapsed 
##   0.240   0.000   0.239 
## set.seed(123); system.time(dat6.1 <- cpda::rplba( n, pVec3.1))
##    user  system elapsed 
##   0.028   0.000   0.028 
## set.seed(123); system.time(dat6.2 <- cpda::rplba( n, pVec3.2))
##    user  system elapsed 
##   0.028   0.000   0.028 
## set.seed(123); system.time(dat6.3 <- cpda::rplba( n, pVec3.3))
##    user  system elapsed 
##   0.028   0.000   0.028

## for ggplot2 and lattice 
tmp5.1 <- data.frame(choice=factor(dat5.1[,1]), rt=dat5.1[,2])  
tmp5.2 <- data.frame(choice=factor(dat5.2[,1]), rt=dat5.2[,2])
tmp5.3 <- data.frame(choice=factor(dat5.3[,1]), rt=dat5.3[,2])
tmp6.1 <- data.frame(choice=factor(dat6.1[,1]), rt=dat6.1[,2])
tmp6.2 <- data.frame(choice=factor(dat6.2[,1]), rt=dat6.2[,2])
tmp6.3 <- data.frame(choice=factor(dat6.3[,1]), rt=dat6.3[,2])


## Below is set-up for ggplot2 and lattice
tmp5.1$fun <- "R" 
tmp5.2$fun <- "R"
tmp5.3$fun <- "R"
tmp6.1$fun <- "C"
tmp6.2$fun <- "C"
tmp6.3$fun <- "C"

tmp5.1$vec <- "1"
tmp5.2$vec <- "2"
tmp5.3$vec <- "3"
tmp6.1$vec <- "1"
tmp6.2$vec <- "2"
tmp6.3$vec <- "3"

df <- rbind(tmp5.1, tmp5.2, tmp5.3, tmp6.1, tmp6.2, tmp6.3)
df$fun <- factor(df$fun)
        
## Show R and C functions produce almost identical distributions        
## Not run (the user needs to install lattice or ggplot2 to plot figures):  
## Set up a colour palette 
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
        "#D55E00", "#CC79A7")
        
require(ggplot2)
p0 <- ggplot(data=df, aes(x = rt, fill=fun, color=fun)) +
  geom_density(alpha=0.2) +
  facet_grid(vec~ choice) +
  scale_fill_manual(values=cb)

## Densities from R function and from C function are overlapped, so
## one is obscured by the other. 
## So even their random engines are not allowed to match numbers exactly
## We know both R and C versions create very similar random numbers
png(filename="r_vs_c.png", width=1024, height=768)
p0
dev.off()

## Lattice and base graphics can plot similar figures, too.   
require(lattice)
histogram( ~rt | vec+choice+fun, data=df, breaks="fd", type="density", 
           xlab="Response Time (s)", 
           panel=function(x, ...) {
                 panel.histogram(x, ...)
                 panel.densityplot(x, darg=list(kernel="gaussian"),...)
  })

## End(Not run)
             
par(mfrow=c(3,2))
hist(tmp5.1[tmp5.1$choice==1,"rt"], breaks="fd", col="gray", freq=FALSE, 
       xlab="RT (s)", main="pLBA-Choice 1")
lines(density(tmp6.1[tmp6.1$choice==1,"rt"]), col="red", lty="dashed",  lwd=1.5)
  
hist(tmp5.1[tmp5.1$choice==2,"rt"], breaks="fd", col="gray", freq=FALSE, 
       xlab="RT (s)", main="pLBA-Choice 2")
lines(density(tmp6.1[tmp6.1$choice==2,"rt"]), col="red", lty="dashed",  lwd=1.5)
  
#######10
hist(tmp5.2[tmp5.2$choice==1,"rt"], breaks="fd", col="gray", freq=FALSE, 
       xlab="RT (s)", main="pLBA-Choice 1")
lines(density(tmp6.2[tmp6.2$choice==1,"rt"]), col="red", lty="dashed",  lwd=1.5)
    
hist(tmp5.2[tmp5.2$choice==2,"rt"], breaks="fd", col="gray", freq=FALSE, 
         xlab="RT (s)", main="pLBA-Choice 2")
lines(density(tmp6.2[tmp6.2$choice==2,"rt"]), col="red", lty="dashed",  lwd=1.5)
    
#######10
hist(tmp5.3[tmp5.3$choice==1,"rt"], breaks="fd", col="gray", freq=FALSE, 
         xlab="RT (s)", main="pLBA-Choice 1")
lines(density(tmp6.3[tmp6.3$choice==1,"rt"]), col="red", lty="dashed",  lwd=1.5)
      
hist(tmp5.3[tmp5.3$choice==2,"rt"], breaks="fd", col="gray", freq=FALSE, 
           xlab="RT (s)", main="pLBA-Choice 2")
lines(density(tmp6.3[tmp6.3$choice==2,"rt"]), col="red", lty="dashed",  lwd=1.5)
par(mfrow=c(1,1))

###########################30
## Example 2               30
## Just another example    30
###########################30
pVec1 <- c(A=1.51, b=2.7, v1=3.32, v2=2.24,  w1=1.51,  w2=3.69,  
           sv=1, rD=0.31, swt=0.5, t0=0.08)
 
pVec2 <- c(A1=1.51, A2=1.51, b1=2.7, b2=2.7, v1=3.32, v2=2.24,
           w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.31, 
           swt=0.5, t0=0.08)
           
system.time(dat1 <- cpda::rplba1( n, pVec1))
system.time(dat2 <- cpda::rplba2( n, pVec2))
system.time(dat3 <- cpda::rplbaR1(n, pVec1))
system.time(dat4 <- cpda::rplbaR2(n, pVec2))

tmp1 <- data.frame(choice=factor(dat1[,1]), rt=dat1[,2])
tmp2 <- data.frame(choice=factor(dat2[,1]), rt=dat2[,2])
tmp3 <- data.frame(choice=factor(dat3[,1]), rt=dat3[,2])
tmp4 <- data.frame(choice=factor(dat4[,1]), rt=dat4[,2])
tmp1$fun <- "rplba1"
tmp2$fun <- "rplba2"
tmp3$fun <- "rplba1-R"
tmp4$fun <- "rplba2-R"
tmp0 <- rbind(tmp1, tmp2, tmp3, tmp4)
tmp0$fun <- factor(tmp0$fun)
  
## Not run:    
require(ggplot2)
ggplot(data = tmp0, aes(x = rt, fill=fun, color=fun)) +
    geom_density(alpha=0.2) +
    facet_grid(.~ choice) +
    scale_fill_manual(values=cb)

## End(Not run)
