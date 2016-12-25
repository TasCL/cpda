library(cpda)
n <- 1e5
pVec1 <- c(A=1.51, b=2.7, v1=3.32, v2=2.24,  w1=1.51,  w2=3.69,  
           sv=1, rD=0.31, swt=0.5, t0=0.08)

pVec2 <- c(A1=1.51, A2=1.51, b1=2.7, b2=2.7, v1=3.32, v2=2.24,
           w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.31, 
           swt=0.5, t0=0.08)

pVec3.1 <- c(A1=1.51, A2=1.51, B1=1.2, B2=1.2, C1=.3, C2=.3, v1=3.32, v2=2.24,
           w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.1, tD=.1,
           swt=0.5, t0=0.08)
pVec3.2 <- c(A1=1.51, A2=1.51, B1=1.2, B2=1.2, C1=.3, C2=.3, v1=3.32, v2=2.24,
           w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.1, tD=.15,
           swt=0.5, t0=0.08)
pVec3.3 <- c(A1=1.51, A2=1.51, B1=1.2, B2=1.2, C1=.3, C2=.3, v1=3.32, v2=2.24,
           w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.15, tD=.1,
           swt=0.5, t0=0.08)


set.seed(123)
system.time(dat1 <- cpda::rplba1(n, pVec1))

set.seed(123)
system.time(dat2 <- cpda::rplba2(n, pVec2))

set.seed(123)
system.time(dat3 <- cpda::rplbaR1(n, pVec1))

set.seed(123)
system.time(dat4 <- cpda::rplbaR2(n, pVec2))

set.seed(123)
system.time(dat5 <- cpda::rplbaR3(n, pVec3))
head(dat5)

all(dat1[,1] == dat2[,1])
all(dat1[,2] == dat2[,2])

all(dat1[,1] == dat3[,1])
all(dat1[,2] == dat3[,2])

all(dat2[,1] == dat4[,1])
all(dat2[,2] == dat4[,2])

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


n <- 1e5
set.seed(123); system.time(dat5.1 <- cpda::rplbaR3(n, pVec3.1))
set.seed(123); system.time(dat5.2 <- cpda::rplbaR3(n, pVec3.2))
set.seed(123); system.time(dat5.3 <- cpda::rplbaR3(n, pVec3.3))
set.seed(123); system.time(dat6.1 <- cpda::rplba3(n, pVec3.1))
set.seed(123); system.time(dat6.2 <- cpda::rplba3(n, pVec3.2))
set.seed(123); system.time(dat6.3 <- cpda::rplba3(n, pVec3.3))

tmp5.1 <- data.frame(choice=factor(dat5.1[,1]), rt=dat5.1[,2])
tmp5.2 <- data.frame(choice=factor(dat5.2[,1]), rt=dat5.2[,2])
tmp5.3 <- data.frame(choice=factor(dat5.3[,1]), rt=dat5.3[,2])
tmp6.1 <- data.frame(choice=factor(dat6.1[,1]), rt=dat6.1[,2])
tmp6.2 <- data.frame(choice=factor(dat6.2[,1]), rt=dat6.2[,2])
tmp6.3 <- data.frame(choice=factor(dat6.3[,1]), rt=dat6.3[,2])


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

tmp0 <- rbind(tmp5.1, tmp5.2, tmp5.3, tmp6.1, tmp6.2, tmp6.3)
tmp0$fun <- factor(tmp0$fun)



library(gridExtra)
library(lattice)
histogram( ~rt | fun + choice , data = tmp0, breaks="fd", type="density", 
                 xlab="Response Time (s)", 
                 panel=function(x, ...) {
                     panel.histogram(x, ...)
                     panel.densityplot(x, darg=list(kernel="gaussian"),...)
                   })


cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
        "#0072B2", "#D55E00", "#CC79A7")
library(ggplot2)
ggplot(data = tmp0, aes(x = rt, fill=fun, color=fun)) +
  geom_density(alpha=0.2) +
  facet_grid(vec~ choice) +
  scale_fill_manual(values=cb)

# hist(tmp1$rt, breaks="fd", freq=FALSE)
# lines(density(tmp1$rt), col="red", lwd=2)

