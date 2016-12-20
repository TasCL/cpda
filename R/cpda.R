#' Probability Density Approximation
#'
#' This package implements probability density approximation in C++ and Open MP
#'
#' @keywords cpda
#' @name cpda
#' @docType package
#' @author  Yi-Shin Lin <yishin.lin@utas.edu.au> \cr
#' Andrew Heathcote <andrew.heathcote@utas.edu.au> \cr
#' William Holmes <william.holmes@vanderbilt.edu>
#' @references Holmes, W. (2015). A practical guide to the Probability Density
#' Approximation (PDA) with improved implementation and error characterization.
#' \emph{Journal of Mathematical Psychology}, \bold{68-69}, 13--24,
#' doi: http://dx.doi.org/10.1016/j.jmp.2015.08.006.
#' @importFrom Rcpp evalCpp
#' @useDynLib cpda
NULL

#' Point-wise Probability Density Approximation
#'
#' \code{logLik_pw} takes individual data point and feeds them sequentially
#' (or in parallel via Open MP) into gaussian kernel to get probability 
#' density directly. 
#'
#' @param y a vector storing empirical data (e.g., RTs)
#' @param yhat a vector storing simulated data (e.g., simualted RTs,
#' using a LBA model).
#' @param h kernel bandwidth. Default value is 0.8 times Silverman's Rule of
#' Thumb, based on simulated data (i.e., yhat). Otherwise, a user entered 
#' value (multiplied by m) will be used. 
#' @param m a bandwidth multiplier. Default is 0.8. Enter 1, uf you do not 
#' @param parallel a switch for parallel processing. Default is FALSE.
#' @return a vector storing log-likelihood for each data point 
#' @examples 
#' #########################
#' ## Example 1           ##
#' #########################
#' ## Test whether logLik_pw match dnorm(0, log=TRUE)
#' rm(list=ls())
#' cpda::logLik_pw(0, rnorm(1e6)) ## -0.9228587
#' dnorm(0, log=TRUE)             ## -0.9189385
#' 
#' h <- 0.8*bw.nrd0(rnorm(1e6));            ## h==0.04541785          
#' cpda::logLik_pw(0, rnorm(1e6), h=h)      ## -0.9196711
#' cpda::logLik_pw(0, rnorm(1e6), h=h, m=1) ## -0.9181807
#' 
#' #########################
#' ## Example 2           ##
#' #########################
#' ## Demo how to use logLik_pw to get pLBA likelihoods 
#' library(cpda); data(lba)
#' str(plba)
#' ## List of 4
#' ## $ DT1 : num [1:695] 0.488 0.801 0.376 0.507 0.532 ...
#' ## $ DT2 : num [1:305] 0.538 0.77 0.568 0.271 0.881 ...
#' ## $ eDT1: num [1:7020] 0.475 0.346 0.42 0.401 0.368 ...
#' ## $ eDT2: num [1:2980] 0.703 0.693 0.704 0.462 0.468 ...
#' 
#' ## Use pointwise pda to get likelihoods for each data point
#' ## This algorithm calculates gaussian kernel directly
#' ## (1) First argument, plba$DT1, is the data.
#' ## (2) Second argument, plba$eDT1, is the simulation.
#' ## (3) The outputs are log-likelihoods. 
#' 
#' outputs <- logLik_pw(plba$DT1, plba$eDT1)
#' sum(outputs) ## Get summed, logged likelihood, 278.6095
#' 
#' #########################
#' ## Example 3           ##
#' #########################
#' rm(list=ls())
#' x <- seq(-3, 3, length.out=100) ## Data
#' samp <- rnorm(1e6)
#' h <- 0.8*stats::bw.nrd0(samp); h          ## Define bin
#' system.time(pw1  <- logLik_pw(x, samp))
#' system.time(pw2  <- logLik_pw(x, samp, h=h, m=0.8))
#' system.time(pw3  <- logLik_pw(x, samp, h=h, m=1))
#' 
#' plot(x, pw1,type="l", lty="dotted")
#' lines(x, pw2, col="darkgreen", lty="dashed")
#' lines(x, pw3, col="blue", lty="dotdash")
#' lines(x, dnorm(x, log=TRUE), col="red")
#' 
#' @export
logLik_pw <- function(y, yhat, h=NA, m=0.8, parallel=FALSE) {
  ## logLik_pw(0, rnorm(1e6)) ## error
  y    <- as.double(sort(y))
  yhat <- as.double(sort(yhat))
  ny   <- as.integer(length(y))
  ns   <- as.integer(length(yhat))
  if(is.na(h)) {
    h_ <- as.double(m*stats::bw.nrd0(yhat))  
  } else {
    h_ <- as.double(m*h)
  }
  
  if (parallel) {
    answer <- .C("logLik_pw_omp", y, yhat, ny, ns, h_,
                 numeric(length(y)), PACKAGE='cpda')
  } else {
    answer <- .C("logLik_pw", y, yhat, ny, ns, h_,
                 numeric(length(y)), PACKAGE='cpda')
  }
  
  return(answer[[6]])
}


