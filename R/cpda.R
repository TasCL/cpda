#' Probability Density Approximation
#'
#' This package implements probability density approximation in C++ 
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

# logLik_pw <- function(y, yhat, h=NA, m=0.8, n=0, parallel=FALSE) {
#   y    <- as.double(y)
#   yhat <- as.double(yhat)
#   ny   <- as.integer(length(y))
#   ns   <- as.integer(length(yhat))
#   h_   <- ifelse(is.na(h), as.double(0), as.double(m*h))
#   m_   <- as.double(m)
#   n_   <- as.double(n)
# 
#   if (parallel) {
#     out <- .C("logLik_pw_omp", y, yhat, ny, ns, h_, m_, n_,
#                  numeric(length(y)), PACKAGE='cpda')
#   } else {
#     out <- .C("logLik_pw", y, yhat, ny, ns, h_, m_, n_,
#                  numeric(length(y)), PACKAGE='cpda')
#   }
#   
#   # return(out[[8]][y_idx])
#   return(out[[8]])
# }


