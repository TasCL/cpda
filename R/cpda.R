#' Probability Density Approximation
#'
#' Implement C++ enhanced probability density approximation
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
#' Instead of using KDE-FFT, \code{logLik_pw} takes individual data point and 
#' feeds them into gaussian kernel to get probability density directly. 
#'
#' @param y a vector storing empirical data (e.g., RTs)
#' @param yhat a vector storing simulated data (e.g., simualted RTs, using a 
#' LBA model).
#' @param h kernel bandwidth. By default 0.8 times Silverman's Rule of Thumb 
#' Bandwidth
#' @param pointwise a boolean switch to turn off pointwise calculation. Should
#' always keep it TRUE. 
#' @return Log-likelihood
#' @examples 
#' library(cpda); data(lba)
#' ## data.table::data.table(d)
#' str(plba)
#' ## List of 4
#' ## $ DT1 : num [1:695] 0.488 0.801 0.376 0.507 0.532 ...
#' ## $ DT2 : num [1:305] 0.538 0.77 0.568 0.271 0.881 ...
#' ## $ eDT1: num [1:7020] 0.475 0.346 0.42 0.401 0.368 ...
#' ## $ eDT2: num [1:2980] 0.703 0.693 0.704 0.462 0.468 ...
#' 
#' ## Use piecewise pda logLik to get likelihood for each data point
#' ## This algorithm calculates gaussian kernel directly
#' ## (1) First argument, plba$DT1, is the data.
#' ## (2) Second argument, plba$eDT1, is the simulation.
#' ## (3) The outputs are likelihoods (not log yet) for each data point
#' 
#' outputs <- logLik_pw(plba$DT1, plba$eDT1)
#' str(outputs)
#' sum(outputs)
#' sum(log(outputs))  ## 278.6095
#' 
#' @export
logLik_pw <- function(y, yhat, h=NA, pointwise=TRUE) {
  y    <- as.double(sort(y))
  yhat <- as.double(sort(yhat))
  ny   <- as.integer(length(y))
  ns   <- as.integer(length(yhat))
  if(is.na(h)) {
    h_ <- as.double(0.8*stats::bw.nrd0(y))
  } else {
    h_ <- as.double(h)
  }
  pw_ <- as.logical(pointwise)
  nt_ <- as.integer(2)
  answer <- .C("logLik_pw", y, yhat, ny, ns, h_,
                 numeric(length(y)), PACKAGE='cpda')
  return(answer[[6]])
}

