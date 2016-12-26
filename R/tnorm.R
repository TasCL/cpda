#' Truncated Normal Distribution
#'
#' The two wrapper functions call two corresponding Rcpp functions to speed up
#' computation. The codes are based on Christopher Jackson's \pkg{msm}
#' (1.5) package and Jonathan Olmsted's \pkg{RcppTN} (0.1-8) package.
#'
#' \code{dtnorm} calculates probability density for a truncated normal
#' distribution with mean equal to \code{mean} and standard deviation equal to
#' \code{sd} and truncates at \code{lower} and \code{upper} range.
#'
#' \code{rtnorm} generates a random number from a truncated Normal
#' distribution with \code{mean} and \code{sd}.
#'
#' @param x,n x in \code{dtnorm} is a vector of quantiles; n in
#' \code{rtnorm} indicates how many random numbers to generate
#' @param mean a scalar or a vector of means
#' @param sd a scalar or a vector of standard deviations
#' @param lower lower truncation value Default as -Inf (ie regular Gaussian)
#' @param upper upper truncation value. Default as Inf (ie regular Gaussian).
#' @param log default 0 (FALSE). Enter 1 to make it TRUE. This will return
#' logged likelihood. This argumement is used only in \code{dtnorm}.
#' @param checks a switch to turn on check functions to check inputs and ouput.
#' Default to FALSE.
#' @export
#' @examples
#' ## Use dtnorm and rtnorm with their default values
#' dtnorm()
#' rtnorm()
#'
#' ## A similar curve plotting example extracted from dnorm functions
#' plot(function(x) dnorm(x, log = FALSE), -2.5, 2.5,
#'      main = "Normal Distribution", ylim=c(0,0.45), ylab="Density")
#' curve(dtnorm(x, lower=-2, upper=2), add=TRUE, col="tomato", lwd=2)
#' mtext("dnorm(x)", adj = 0)
#' mtext("dtnorm(x)", col = "tomato", adj = 1)
dtnorm <- function(x=0, mean=0, sd=1, lower=-Inf, upper=Inf, log=0, checks=FALSE)
{
  mean  <- rep(mean,  length=length(x))
  sd    <- rep(sd,    length=length(x))
  lower <- rep(lower, length=length(x))
  upper <- rep(upper, length=length(x))
  log   <- rep(log, length=length(x))
  
  if (checks) { checkInputs(mean, sd, lower, upper) }
  
  ## inside C++, I use islog as the variable name for the boolean log, so
  ## I can use cmath's log function
  out <- .Call("dtn_wrapper", x_ = x,
               mean_ = mean,   sd_ = sd,
               lower_ = lower, upper_ = upper,
               islog_=log)
  
  if (checks) { checkOutputs(out) }
  return(out)
}

#' @export
#' @rdname dtnorm
rtnorm <- function(n=1, mean=0, sd=1, lower=-Inf, upper=Inf, checks=FALSE)
{
  ## The size of mean will determine the draw numbers
  if (length(n) > 1) n <- length(n)
  mean  <- rep(mean,  length=n)
  sd    <- rep(sd,    length=n)
  lower <- rep(lower, length=n)
  upper <- rep(upper, length=n)
  
  if (checks) { checkInputs(mean, sd, lower, upper) }
  
  out <- .Call("rtn_wrapper", mean_ = mean,   sd_ = sd,
               lower_ = lower, upper_ = upper)
  
  if (checks) { checkOutputs(out) }
  return(out)
}

checkInputs <- function(m, s, l, u)
{
  if (!(length(m) == length(s) & length(m) == length(l) & length(m) == length(u)))
  {
    stop("Input vectors are not all same length. Nothing done.")
  }
}

checkOutputs <- function(out)
{
  if (any(is.na(out)))
  {
    warning("NAs returned. Check for invalid parameters.")
  }
}
