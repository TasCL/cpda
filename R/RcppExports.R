# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @export
gaussian <- function(x, xtilde, h) {
    .Call('_cpda_gaussian', PACKAGE = 'cpda', x, xtilde, h)
}

#' @export
lik_fft <- function(y, yhat, h = 0, m = 0.8, p = 10, n = 0L) {
    .Call('_cpda_lik_fft', PACKAGE = 'cpda', y, yhat, h, m, p, n)
}

#' Point-wise Probability Density Approximation
#'
#' \code{lik_pw} takes each observation and sequentially (or concurrently)
#' via Open MP) conduct KDEs. This function is for testing purpose. If you 
#' wish to conduct KDE smoothing point-by-point, use 'density'
#' in \code{stats}, which uses a similar FFT algorithm as \code{lik_fft}. 
#'
#' @param y a vector storing empirical observations (e.g., RTs).
#' @param yhat a vector storing simulations.
#' @param h kernel bandwidth. Default value is 0.8 times Silverman's Rule of
#' Thumb, based on simulated data (i.e., yhat).
#' @param m a bandwidth multiplier. Default is 0.8.
#' @param n the number of simulation. When n=0, where the function will count
#' the numbers of observation in the simulated histogram. If simulating 
#' a defective distribution, one should enter the simulation numbers. 
#' @param parallel a switch for parallel processing via OMP. Default is FALSE. 
#' @return a vector storing log-likelihoods
#' @references 
#' Holmes, W. (2015). A practical guide to the Probability Density
#' Approximation (PDA) with improved implementation and error characterization.
#' \emph{Journal of Mathematical Psychology}, vol. 68-69, 13--24,
#' doi: \url{http://dx.doi.org/10.1016/j.jmp.2015.08.006}. \cr
#' 
#' Turner, B. M. & Sederberg, P. B. (2014). A generalized, likelihood-free 
#' method for posterior estimation. \emph{Psychonomic Bulletin Review}, 21(2), 
#' 227-250. doi: \url{http://dx.doi.org/10.3758/s13423-013-0530-0}. \cr
#' 
#' Brown, S. & Heathcote, A. (2008). The simplest complete model of choice 
#' response time: Linear ballistic accumulation. \emph{Cognitive Psychology}, 
#' 57, 153-178. doi: \url{http://dx.doi.org/10.1016/j.cogpsych.2007.12.002}.
#' @export
#' @examples
#' ## Example 1 tests if Lik_pw match 'dnorm' and shows how to adjust 
#' ## bandwidth  
#' rm(list=ls())
#' lik_pw(0, rnorm(1e6)) ## 0.3974386
#' dnorm(0)              ## 0.3989423
#'
#' h <- 0.8*bw.nrd0(rnorm(1e6));            ## 0.04542646 
#' lik_pw(0, rnorm(1e6), h=h)      ## 0.3996198
#' lik_pw(0, rnorm(1e6), h=h, m=1) ## 0.3985692
#'
#' ## Example 2 demostrates how to use Lik_pw to get pLBA likelihoods
#' data(lba)
#' str(plba)
#' ## List of 4
#' ## $ DT1 : num [1:695] 0.488 0.801 0.376 0.507 0.532 ...
#' ## $ DT2 : num [1:305] 0.538 0.77 0.568 0.271 0.881 ...
#' ## $ eDT1: num [1:7020] 0.475 0.346 0.42 0.401 0.368 ...
#' ## $ eDT2: num [1:2980] 0.703 0.693 0.704 0.462 0.468 ...
#'
#' ## Use pointwise pda to get likelihoods for each data point
#' ## This algorithm calculates via a standard gaussian kernel directly
#' ## (1) First argument, plba$DT1, is the data.
#' ## (2) Second argument, plba$eDT1, is the simulation.
#' ## (3) The output is likelihood.
#'
#' output <- lik_pw(plba$DT1, plba$eDT1)
#' sum(output) ## Get summed, logged likelihood
#'
#' #########################
#' ## Example 3           ##
#' #########################
#' rm(list=ls())
#' n      <- 1e5
#' x      <- seq(-3,3, length.out=100) ## Support
#' xlabel <- "Observations" 
#' ylabel <- "Density" 
#' 
#' ## Approximate Gaussian densities
#' samp <- rnorm(n)
#' pw1  <- lik_pw(x, samp)
#' pw2  <- approx(density(samp)$x, density(samp)$y, x)$y
#' plot(x,  pw1, type="l", lty="longdash", xlab=xlabel, ylab=ylabel)
#' lines(x, pw2, lwd=1.5, lty="dotted")
#' lines(x, dnorm(x), lwd=2)
#' 
#' samp <- gamlss.dist::rexGAUS(n, mu=-2, sigma=1, nu=1)
#' system.time(pw1 <- lik_pw(x, samp, parallel=TRUE))
#' system.time(pw2 <- approx(density(samp)$x, density(samp)$y, x)$y)
#' plot(x,  pw1, type="l", lty="longdash", xlab=xlabel, ylab=ylabel)
#' lines(x, pw2, lwd=1.5, lty="dotted")
#' lines(x, gamlss.dist::dexGAUS(x, mu=-2, sigma=1, nu=1), col="red")
#'
#' ## Approximate densities of linear regression with Gaussian noise: 
#' ## y = ax + b + N(0, s)
#' theta <- c(a=7.5, b=3.5, s=5)
#' y     <- rnorm(length(x), mean=theta[2]+theta[1]*x, sd=theta[3])
#' dat   <- cbind(x, y)
#' plot(x, y, main="Linear Regression")
#' 
#' ## -- Because means for each data point differ, we need to generate 
#' ## simulation (ie samp) for each data point. That is, the likelihood for  
#' ## each data point is calculated by different simulated likelihood functions
#' ## -- We use sapply to gain a little speedy up. 
#' ## -- 'density' function cannot calculate density for only 1 data point 
#' pw1 <- sapply(x, function(s) {
#'   samp  <- rnorm(n, mean=theta[2]+theta[1]*s, sd=theta[3])
#'   lik_pw(s, samp)
#' })
#' plot(x,  pw1,  type="l", lty="longdash", xlab=xlabel, ylab=ylabel)
#' lines(x,  dnorm(x, theta[2]+theta[1]*x, theta[3]), col="red")
#' 
#' #########################
#' ## Example 4: LBA      ##
#' #########################
#' rm(list=ls())
#' ## Assuming that this is empirical data
#' y    <- rtdists::rLBA(1e3, A=.5, b=1, t0=.25, mean_v=c(2.4, 1.2), sd_v=c(1,1.2))
#' rt1  <- y[y$response==1,"rt"]
#' rt2  <- y[y$response==2,"rt"]
#' srt1 <- sort(rt1)
#' srt2 <- sort(rt2)
#' summary(rt1); summary(rt2)
#'     
#' n <- 1e5
#' pvec <- c(b=1, A=.5, v1=2.4, v2=1.2, sv1=1, sv2=1.2, t0=.25)
#' samp <- cpda::rlba1(n, pvec)
#' samp1 <- samp[samp[,2]==1, 1]
#' samp2 <- samp[samp[,2]==2, 1]
#' pw1 <- lik_pw(srt1, samp1, n=n)
#' pw2 <- lik_pw(srt2, samp2, n=n)
#'     
#' den0 <- rtdists::dLBA(y$rt, y$response, A=.5, b=1, t0=.25, 
#'   mean_v=c(2.4, 1.2), sd_v=c(1, 1.2))
#' df0 <- cbind(y, den0)
#' df1 <- df0[df0[,2]==1,]
#' df2 <- df0[df0[,2]==2,]
#' den1 <- df1[order(df1[,1]),3]
#' den2 <- df2[order(df2[,1]),3]
#'   
#' plot(srt1,  den1, type="l")
#' lines(srt2, den2)
#' lines(srt1, pw1, col="red")
#' lines(srt2, pw2, col="red")
#'     
#' #########################
#' ## Example 5           ##
#' #########################
#' n   <- 1e5
#' x   <- seq(-3,3, length.out=100) ## Support
#' sam <- rnorm(n)
#' 
#' ## Tested on 12-core CPU
#' rbenchmark::benchmark(replications=rep(10, 3),
#'    pw       <- cpda::lik_pw(x, sam),
#'    pw_omp   <- cpda::lik_pw(x, sam, parallel = T),
#'    columns=c('test', 'elapsed', 'replications'))
#' ##                                           test elapsed replications
#' ## 1                   pw <- cpda::lik_pw(x, sam)   2.570           10
#' ## 3                   pw <- cpda::lik_pw(x, sam)   2.485           10
#' ## 5                   pw <- cpda::lik_pw(x, sam)   2.484           10
#' ## 2 pw_omp <- cpda::lik_pw(x, sam, parallel = T)   0.993           10
#' ## 4 pw_omp <- cpda::lik_pw(x, sam, parallel = T)   1.119           10
#' ## 6 pw_omp <- cpda::lik_pw(x, sam, parallel = T)   1.024           10
#' 
#' @export
lik_pw <- function(x, xtilde, h = 0, m = 0.8, n = 0L) {
    .Call('_cpda_lik_pw', PACKAGE = 'cpda', x, xtilde, h, m, n)
}

#' @export
n1PDF <- function(x, nsim, b, A, mean_v, sd_v, t0, h_in = 0, k = 0.09, debug = FALSE) {
    .Call('_cpda_n1PDF', PACKAGE = 'cpda', x, nsim, b, A, mean_v, sd_v, t0, h_in, k, debug)
}

rtn_scalar <- function(mean, sd, l, u) {
    .Call('_cpda_rtn_scalar', PACKAGE = 'cpda', mean, sd, l, u)
}

#' @export
dtn_ <- function(x, mean, sd, lower, upper, lp) {
    .Call('_cpda_dtn_', PACKAGE = 'cpda', x, mean, sd, lower, upper, lp)
}

#' @export
dtn <- function(x, mean, sd, lower, upper, log = FALSE) {
    .Call('_cpda_dtn', PACKAGE = 'cpda', x, mean, sd, lower, upper, log)
}

#' @export
ptn_ <- function(q, mean, sd, lower, upper, lt, lp) {
    .Call('_cpda_ptn_', PACKAGE = 'cpda', q, mean, sd, lower, upper, lt, lp)
}

#' @export
ptn <- function(q, mean, sd, lower, upper, lt = TRUE, lp = FALSE) {
    .Call('_cpda_ptn', PACKAGE = 'cpda', q, mean, sd, lower, upper, lt, lp)
}

#' @rdname dtnorm
#' @export
rtnorm <- function(n, mean, sd, lower, upper) {
    .Call('_cpda_rtnorm', PACKAGE = 'cpda', n, mean, sd, lower, upper)
}

#' @export
dexGAUS <- function(x, mu, sigma, tau, log = FALSE) {
    .Call('_cpda_dexGAUS', PACKAGE = 'cpda', x, mu, sigma, tau, log)
}

#' @export
pexGAUS <- function(q, mu, sigma, tau, lt = TRUE, lp = FALSE) {
    .Call('_cpda_pexGAUS', PACKAGE = 'cpda', q, mu, sigma, tau, lt, lp)
}

#' @export
qexGAUS <- function(p, mu, sigma, tau, lt = TRUE, lp = FALSE) {
    .Call('_cpda_qexGAUS', PACKAGE = 'cpda', p, mu, sigma, tau, lt, lp)
}

#' @export
dinvGAUS <- function(x, mu, lambda, log = FALSE) {
    .Call('_cpda_dinvGAUS', PACKAGE = 'cpda', x, mu, lambda, log)
}

#' @export
pinvGAUS <- function(q, mu, lambda, lt = TRUE, lp = FALSE) {
    .Call('_cpda_pinvGAUS', PACKAGE = 'cpda', q, mu, lambda, lt, lp)
}

#' @export
qinvGAUS <- function(p, mu, lambda, lt = TRUE, lp = FALSE) {
    .Call('_cpda_qinvGAUS', PACKAGE = 'cpda', p, mu, lambda, lt, lp)
}

#' @export
rinvGAUS_ <- function(n, mu, lambda) {
    .Call('_cpda_rinvGAUS_', PACKAGE = 'cpda', n, mu, lambda)
}

#' @export
rinvGAUS <- function(n, mu, lambda) {
    .Call('_cpda_rinvGAUS', PACKAGE = 'cpda', n, mu, lambda)
}

#' @export
InvG <- function(x) {
    invisible(.Call('_cpda_InvG', PACKAGE = 'cpda', x))
}

#' @export
rlba <- function(n, b, A, mean_v, sd_v, t0) {
    .Call('_cpda_rlba', PACKAGE = 'cpda', n, b, A, mean_v, sd_v, t0)
}

#' @export
rlba_n1 <- function(n, b, A, mean_v, sd_v, t0) {
    .Call('_cpda_rlba_n1', PACKAGE = 'cpda', n, b, A, mean_v, sd_v, t0)
}

#' Generate Random Choice Response Times using pLBA Model
#'
#' This function uses two-accumulator piecewise LBA model to generate random
#' choice RTs. There are 3 variants: \code{rplba}, \code{rplba1}, and
#' \code{rplba2}. Each of them has a corresponding R version,
#' \code{rplbaR}, \code{rplbaR1}, and \code{rplbaR2}, for the purpose of
#' speed testing. Because the difference of random number generators in C and
#' R, they do not generate exactly identical RTs.  When generating large
#' enough observations, the distributions generated by R and C will match.
#'
#' The main function \code{rplba} implements a more flexible
#' version of pLBA random number generator than the other two. It uses the
#' following parameterisation (order matter):
#'
#' \itemize{
#' \item \bold{\emph{A1}} accumulator 1 start-point upper bound. \code{A} is
#' the upper bound of the interval \code{[0, A]}, which is used by an uniform
#' distribution to generate a start-point. Average amount of
#' prior evidence (i.e., before accumulation process even begins) across trials
#' is \code{A/2}.
#' \item \bold{\emph{A2}} accumulator 2 start-point upper bound.
#' \item \bold{\emph{B1}} accumulator 1 traveling distance. Note this is not
#' a decision threshold!. LBA convention denotes decision threshold/caution as
#' b (lowercase) and traveling distance as B (uppercase). \code{B=b-A} is
#' the traveling distance, and \code{b-A/2} is a measure of average
#' \emph{decision caution}.
#' \item \bold{\emph{B2}} accumulator 2 traveling distance.
#' \item \bold{\emph{C1}} the amount of traveling distance change for
#' accumulator 1 at the stage 2.
#' \item \bold{\emph{C2}} the amount of traveling distance change for
#' accumulator 2 at the stage 2.
#' \item \bold{\emph{v1}} accumulator 1 drift rate, stage 1
#' \item \bold{\emph{v2}} accumulator 2 drift rate, stage 1
#' \item \bold{\emph{w1}} accumulator 1 drift rate, stage 2
#' \item \bold{\emph{w2}} accumulator 2 drift rate, stage 2
#' \item \bold{\emph{sv1}} accumulator 1 drift rate standard deviation,
#' stage 1.
#' \item \bold{\emph{sv2}} accumulator 2 drift rate standard deviation,
#' stage 1.
#' \item \bold{\emph{sw1}} accumulator 1 drift rate standard deviation,
#' stage 2.
#' \item \bold{\emph{sw2}} accumulator 2 drift rate standard deviation,
#' stage 2.
#' \item \bold{\emph{rD}} the delay duration while stage 1 drift rate switches
#' to stage 2 drift rate
#' \item \bold{\emph{tD}} the delay duration while stage 1 threshold switches
#' to stage 2 threshold
#' \item \bold{\emph{swt}} switch time, usually determined by experimental
#' design.
#' \item \bold{\emph{t0}} non-decision time in second.
#' }
#'
#' \code{rplba1} uses the following parameterisation:
#'
#' \itemize{
#' \item \bold{\emph{A}} a common start-point interval for both accumulators.
#' \item \bold{\emph{b}} a common response threshold for both accumulators.
#' \item \bold{\emph{v1}} accumulator 1 drift rate, stage 1
#' \item \bold{\emph{v2}} accumulator 2 drift rate, stage 1
#' \item \bold{\emph{w1}} accumulator 1 drift rate, stage 2
#' \item \bold{\emph{w2}} accumulator 2 drift rate, stage 2
#' \item \bold{\emph{sv}} a common standard deviation for both accumulators
#' \item \bold{\emph{rD}} a delay period while drift rate switch to a
#' second stage process
#' \item \bold{\emph{swt}} switch time, usually determined by experimental
#' design
#' \item \bold{\emph{t0}} non-decision time in second.
#' }
#'
#' \code{rplba2} uses the following parameterisation:
#'
#' \itemize{
#' \item \bold{\emph{A1}} start-point interval of the accumulator 1.
#' \item \bold{\emph{A2}} start-point interval of the accumulator 2.
#' \item \bold{\emph{b1}} accumulator 1 response threshold.
#' \item \bold{\emph{b2}} accumulator 2 response threshold.
#' \item \bold{\emph{v1}} accumulator 1 drift rate, stage 1
#' \item \bold{\emph{v2}} accumulator 2 drift rate, stage 1
#' \item \bold{\emph{w1}} accumulator 1 drift rate, stage 2
#' \item \bold{\emph{w2}} accumulator 2 drift rate, stage 2
#' \item \bold{\emph{sv1}} the standard deviation of accumulator 1 drirt rate
#' during stage 1.
#' \item \bold{\emph{sv2}} the standard deviation of accumulator 2 drirt rate
#' during stage 1.
#' \item \bold{\emph{sw1}} the standard deviation of accumulator 1 drirt rate
#' during stage 2.
#' \item \bold{\emph{sw2}} the standard deviation of accumulator 2 drirt rate
#' during stage 2.
#' \item \bold{\emph{rD}} a delay period while drift rate switch to a
#' second stage process
#' \item \bold{\emph{swt}} switch time, usually determined by experimental
#' design
#' \item \bold{\emph{t0}} non-decision time in second.
#' }
#'
#' @param n number of observations. Must be an integer
#' @param pVec a numeric vector storing pLBA model parameters. The sequence is
#' critical. See details for the sequence.
#' @return A \code{n x 2} matrix with a first column storing choices and second
#' column storing response times.
#' @examples
#' ################
#' ## Example 1  ##
#' ################
#' pVec3.1 <- c(A1=1.51, A2=1.51, B1=1.2, B2=1.2, C1=.3, C2=.3, v1=3.32,
#'              v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.1,
#'              tD=.1, swt=0.5, t0=0.08)
#' pVec3.2 <- c(A1=1.51, A2=1.51, B1=1.2, B2=1.2, C1=.3, C2=.3, v1=3.32,
#'              v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.1,
#'              tD=.15, swt=0.5, t0=0.08)
#' pVec3.3 <- c(A1=1.51, A2=1.51, B1=1.2, B2=1.2, C1=.3, C2=.3, v1=3.32,
#'              v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.15,
#'              tD=.1, swt=0.5, t0=0.08)
#'
#' n <- 1e5
#' set.seed(123); system.time(dat5.1 <- cpda::rplbaR(n, pVec3.1))
#' set.seed(123); system.time(dat5.2 <- cpda::rplbaR(n, pVec3.2))
#' set.seed(123); system.time(dat5.3 <- cpda::rplbaR(n, pVec3.3))
#' set.seed(123); system.time(dat6.1 <- cpda::rplba( n, pVec3.1))
#' set.seed(123); system.time(dat6.2 <- cpda::rplba( n, pVec3.2))
#' set.seed(123); system.time(dat6.3 <- cpda::rplba( n, pVec3.3))
#' tmp5.1 <- data.frame(choice=factor(dat5.1[,1]), rt=dat5.1[,2])
#' tmp5.2 <- data.frame(choice=factor(dat5.2[,1]), rt=dat5.2[,2])
#' tmp5.3 <- data.frame(choice=factor(dat5.3[,1]), rt=dat5.3[,2])
#' tmp6.1 <- data.frame(choice=factor(dat6.1[,1]), rt=dat6.1[,2])
#' tmp6.2 <- data.frame(choice=factor(dat6.2[,1]), rt=dat6.2[,2])
#' tmp6.3 <- data.frame(choice=factor(dat6.3[,1]), rt=dat6.3[,2])
#'
#' tmp5.1$fun <- "R"
#' tmp5.2$fun <- "R"
#' tmp5.3$fun <- "R"
#' tmp6.1$fun <- "C"
#' tmp6.2$fun <- "C"
#' tmp6.3$fun <- "C"
#'
#' tmp5.1$vec <- "1"
#' tmp5.2$vec <- "2"
#' tmp5.3$vec <- "3"
#' tmp6.1$vec <- "1"
#' tmp6.2$vec <- "2"
#' tmp6.3$vec <- "3"
#'
#' df <- rbind(tmp5.1, tmp5.2, tmp5.3, tmp6.1, tmp6.2, tmp6.3)
#' df$fun <- factor(df$fun)
#'
#' ## Show R and C functions produce almost identical distributions
#' \dontrun{
#' ## Set up a colour palette
#' cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
#'         "#D55E00", "#CC79A7")
#'
#' require(ggplot2)
#' ggplot(data=df, aes(x = rt, fill=fun, color=fun)) +
#'   geom_density(alpha=0.2) +
#'   facet_grid(vec~ choice) +
#'   scale_fill_manual(values=cb)
#'
#' ## Or you can use lattice or base graphics
#' require(lattice)
#' histogram( ~rt | vec+choice+fun, data=df, breaks="fd", type="density",
#'            xlab="Response Time (s)",
#'            panel=function(x, ...) {
#'                  panel.histogram(x, ...)
#'                  panel.densityplot(x, darg=list(kernel="gaussian"),...)
#'   })
#' }
#'
#' par(mfrow=c(3,2))
#' hist(tmp5.1[tmp5.1$choice==1,"rt"], breaks="fd", col="gray", freq=FALSE,
#'        xlab="RT (s)", main="pLBA-Choice 1")
#' lines(density(tmp6.1[tmp6.1$choice==1,"rt"]), col="red", lty="dashed",  lwd=1.5)
#'
#' hist(tmp5.1[tmp5.1$choice==2,"rt"], breaks="fd", col="gray", freq=FALSE,
#'        xlab="RT (s)", main="pLBA-Choice 2")
#' lines(density(tmp6.1[tmp6.1$choice==2,"rt"]), col="red", lty="dashed",  lwd=1.5)
#'
#' #############
#' hist(tmp5.2[tmp5.2$choice==1,"rt"], breaks="fd", col="gray", freq=FALSE,
#'        xlab="RT (s)", main="pLBA-Choice 1")
#' lines(density(tmp6.2[tmp6.2$choice==1,"rt"]), col="red", lty="dashed",  lwd=1.5)
#'
#' hist(tmp5.2[tmp5.2$choice==2,"rt"], breaks="fd", col="gray", freq=FALSE,
#'          xlab="RT (s)", main="pLBA-Choice 2")
#' lines(density(tmp6.2[tmp6.2$choice==2,"rt"]), col="red", lty="dashed",  lwd=1.5)
#'
#' #############
#' hist(tmp5.3[tmp5.3$choice==1,"rt"], breaks="fd", col="gray", freq=FALSE,
#'          xlab="RT (s)", main="pLBA-Choice 1")
#' lines(density(tmp6.3[tmp6.3$choice==1,"rt"]), col="red", lty="dashed",  lwd=1.5)
#'
#' hist(tmp5.3[tmp5.3$choice==2,"rt"], breaks="fd", col="gray", freq=FALSE,
#'            xlab="RT (s)", main="pLBA-Choice 2")
#' lines(density(tmp6.3[tmp6.3$choice==2,"rt"]), col="red", lty="dashed",  lwd=1.5)
#' par(mfrow=c(1,1))
#'
#' ################
#' ## Example 2  ##
#' ################
#' pVec1 <- c(A=1.51, b=2.7, v1=3.32, v2=2.24,  w1=1.51,  w2=3.69,
#'            sv=1, rD=0.31, swt=0.5, t0=0.08)
#'
#' pVec2 <- c(A1=1.51, A2=1.51, b1=2.7, b2=2.7, v1=3.32, v2=2.24,
#'            w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.31,
#'            swt=0.5, t0=0.08)
#'
#' system.time(dat1 <- cpda::rplba1( n, pVec1))
#' system.time(dat2 <- cpda::rplba2( n, pVec2))
#' system.time(dat3 <- cpda::rplbaR1(n, pVec1))
#' system.time(dat4 <- cpda::rplbaR2(n, pVec2))
#'
#' tmp1 <- data.frame(choice=factor(dat1[,1]), rt=dat1[,2])
#' tmp2 <- data.frame(choice=factor(dat2[,1]), rt=dat2[,2])
#' tmp3 <- data.frame(choice=factor(dat3[,1]), rt=dat3[,2])
#' tmp4 <- data.frame(choice=factor(dat4[,1]), rt=dat4[,2])
#' tmp1$fun <- "rplba1"
#' tmp2$fun <- "rplba2"
#' tmp3$fun <- "rplba1-R"
#' tmp4$fun <- "rplba2-R"
#' tmp0 <- rbind(tmp1, tmp2, tmp3, tmp4)
#' tmp0$fun <- factor(tmp0$fun)
#'
#' \dontrun{
#' require(ggplot2)
#' ggplot(data = tmp0, aes(x = rt, fill=fun, color=fun)) +
#'     geom_density(alpha=0.2) +
#'     facet_grid(.~ choice) +
#'     scale_fill_manual(values=cb)
#' }
#' @export
rplba <- function(n, pVec) {
    .Call('_cpda_rplba', PACKAGE = 'cpda', n, pVec)
}

#' @rdname rplba
#' @export
rplba1 <- function(n, pVec) {
    .Call('_cpda_rplba1', PACKAGE = 'cpda', n, pVec)
}

#' @rdname rplba
#' @export
rplba2 <- function(n, pVec) {
    .Call('_cpda_rplba2', PACKAGE = 'cpda', n, pVec)
}

#' @rdname rplba
#' @export
rplba3 <- function(n, pVec) {
    .Call('_cpda_rplba3', PACKAGE = 'cpda', n, pVec)
}

#' @rdname rplba
#' @export
rplba4 <- function(n, pVec) {
    .Call('_cpda_rplba4', PACKAGE = 'cpda', n, pVec)
}

#' @export
rplba5 <- function(n, pVec) {
    .Call('_cpda_rplba5', PACKAGE = 'cpda', n, pVec)
}

#' @export
rplba6 <- function(n, pVec) {
    .Call('_cpda_rplba6', PACKAGE = 'cpda', n, pVec)
}

#' A simple and fast quantile calculator
#'
#' A C++ quantile function.
#'
#' @param y a data vector
#' @param q nth quantile. Enter proportion, such as .25 or .75.
#' @examples
#' y <- rnorm(100)
#' q <- cquantile(y, .25)
#' @export
cquantile <- function(y, q) {
    .Call('_cpda_cquantile', PACKAGE = 'cpda', y, q)
}

#' Silverman's Rule of Thumb Bandwidth for Kernel Density Estimation
#'
#' A C++ version of Silverman's rule-of-thumb bandwidth. This is similar
#' with R's \code{bw.nrd0(x)}
#'
#' @param y a data vector
#' @param m a multiplier to adjust the SROT proportionally.
#'
#' @seealso
#' \code{\link{bw.nrd}}, \code{\link{bandwidth.nrd}}.
#' @export
#' @examples
#' data(lba)
#' h <-cpda::bwNRD0(plba$DT1, 0.8)
#'
bwNRD0 <- function(y, m) {
    .Call('_cpda_bwNRD0', PACKAGE = 'cpda', y, m)
}

#' @export
histc <- function(x, edge) {
    .Call('_cpda_histc', PACKAGE = 'cpda', x, edge)
}

#' @export
histd <- function(yhat, z, n) {
    .Call('_cpda_histd', PACKAGE = 'cpda', yhat, z, n)
}

#' Calculate Histogram Edges 
#'
#' This is an internal function. The user may not use it.
#' 
#' @param z a grid a scalar, usually created by 
#' \code{z = arma::linspace<arma::vec>(z0, z1, 1<<(int)p);}
#' 
#' @examples
#' set.seed(123)
#' dat1 <- stats::rnorm(1e2)
#' h  <- bw.nrd0(dat1)
#' z0 <- min(dat1) - 3*h 
#' z1 <- max(dat1) + 3*h
#' ngrid <- 2^10
#' dt <- (z1 - z0) / (ngrid - 1)
#' z  <- z0 + (0:(ngrid-1)) * dt
#' ## Same as using seq function
#' ## seq(z0, z1, length.out=ngrid)
#' 
#' binedge1 <- c(z - 0.5*dt, z[ngrid] + 0.5*dt)
#' binedge2 <- as.vector(cpda::getEdges(z))
#' all.equal(binedge1, binedge2)
#' mean(binedge - as.vector(tmp))
#' ## [1] 2.640334e-17
#' 
#' ## 4.00 vs. 14.1 microseconds
#' ## library(microbenchmark)
#' ## res <- microbenchmark(
#' ##    binedge1 <- c(z - 0.5*dt, z[ngrid] + 0.5*dt),
#' ##    binedge2 <- as.vector(cpda::getEdges(z)),
#' ##   times=10L)
#' @export
getEdges <- function(z) {
    .Call('_cpda_getEdges', PACKAGE = 'cpda', z)
}

#' @export
getFilter <- function(m, M, h, p) {
    .Call('_cpda_getFilter', PACKAGE = 'cpda', m, M, h, p)
}

