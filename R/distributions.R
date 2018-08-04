#' The ex-Gaussian Distribution 
#' 
#' Density, distribution function, quantile function and random generation for 
#' the ex-Gaussian distribution with \eqn{mean = mu + tau} and 
#' \eqn{standard deviation} \eqn{ = sigma^{2} + tau^{2}}. See Table 1 in 
#' Cousineau, Brown, & Heathcote (2004) for the equation.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu vector of mu parameters. 
#' @param sigma vector of sigma parameters.  
#' @param tau vector of tau parameters.  
#' @param log,log.p  logical; if TRUE, probabilities p are given as \code{log(p)}
#' @param lower.tail logical; if TRUE, (default), probabilities are 
#' \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' 
#' @return \code{dexG} gives the density. \code{pexG} gives the distribtuion
#' function, \code{qexG} gives the quantile function, and \code{rexG} generates
#' random deviates. \eqn{sigma < 0} or \eqn{tau < 0} is an error. 
#' @examples
#' x <- seq(100, 600, length.out = 1024)
#' den <- cpda::dexG(x, mu=300, sigma=35, tau=100)
#' 
#' ## Use Table 3 in Cousineau, Brown, & Heathcote (2004) as an example
#' n <- 1024
#' mu <- 910.56
#' sigma <- 44.721
#' tau <- 89.443
#' sam <- cpda::rexG(n, mu, sigma, tau)
#' hist(sam, freq = FALSE, breaks = "FD", main = "ex-Gaussian Distribution", 
#'   xlab = "RT (ms)", ylab="Density")
#'   
#' @export
dexG <- function(x, mu = 5, sigma = 1, tau = 1, log = FALSE) {
  dexGAUS(x, mu, sigma, tau, log)[,1]
}

#' @rdname dexG
#' @export
pexG <- function(q, mu = 5, sigma = 1, tau = 1, lower.tail = TRUE, 
  log.p = FALSE) {
  pexGAUS(q, mu, sigma, tau, lower.tail, log.p)[,1]
}

#' @rdname dexG
#' @export
qexG <- function(p, mu = 5, sigma = 1, tau = 1, lower.tail = TRUE, 
  log.p = FALSE) {
  qexGAUS(p, mu, sigma, tau, lower.tail, log.p)[,1]
}

#' @rdname dexG
#' @export
rexG <- function(n, mu = 5, sigma = 1, tau = 1) { 
  if (any(sigma <= 0)) stop("sigma must be positive\n") 
  if (any(tau <= 0)) stop("nu must be positive\n") 
  if (any(n <= 0)) stop("n must be a positive integer\n")  
  n <- ceiling(n)
  p <- runif(n)
  return(qexG(p, mu=mu, sigma=sigma, tau=tau))
}

#' The inverse-Gaussian Distribution 
#' 
#' Density, distribution function, quantile function and random generation for 
#' the inverse-Gaussian distribution with \eqn{mean = mu} and 
#' \eqn{standard deviation} \eqn{ = mu^{3} / lambda}. See Table 1 in 
#' Cousineau, Brown, & Heathcote (2004) for the equation. kappa assumes 0.
#' 
#' The inverse-Gaussian is also known as the Wald distribution.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu vector of mu parameters. 
#' @param lambda vector of lambda parameters.  
#' @param log,log.p  logical; if TRUE, probabilities p are given as \code{log(p)}
#' @param lower.tail logical; if TRUE, (default), probabilities are 
#' \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' 
#' @return \code{dinvG} gives the density. \code{pinvG} gives the distribtuion
#' function, \code{qinvG} gives the quantile function, and \code{rinvG} 
#' generates random deviates. \eqn{mu <= 0} or \eqn{lambda <= 0} is an error. 
#' @examples
#' x <- seq(0, 1, length.out = 1024)
#' den <- cpda::dinvG(x, mu=.275, lambda=2)
#' plot(x, den, type = "l", lwd = 2)
#' 
#' ## Use Table 3 in Cousineau, Brown, & Heathcote (2004) as an example
#' n <- 4096
#' sam <- cpda::rinvG(n, 275, 2000) + 725
#' hist(sam, freq = FALSE, breaks = "FD", main = "Wald Distribution", 
#'   xlab = "RT (ms)", ylab="Density")
#'   
#' @export
dinvG <- function(x, mu = 1, lambda = 1, log = FALSE) {
  dinvGAUS(x, mu, lambda, log)[,1]
}

#' @rdname dinvG
#' @export
pinvG <- function(q, mu = 1, lambda = 1, lower.tail = TRUE, 
  log.p = FALSE) {
  pinvGAUS(q, mu, lambda, lower.tail, log.p)[,1]
}

#' @rdname dinvG
#' @export
qinvG <- function(p, mu = 1, lambda = 1, lower.tail = TRUE, log.p = FALSE) {
  qinvGAUS(p, mu, lambda, lower.tail, log.p)[,1]
}

#' @rdname dinvG
#' @export
rinvG <- function(n, mu = 1, lambda = 1) { 
  if (length(n) != 1) stop("n must be a scalar\n")
  if (n <= 0) stop("n must be a positive integer\n")  
  if (length(mu) != length(lambda)) stop ("mu and lambda must have same number of elements")

  if (length(mu) == 1) {
    out <- rinvGAUS_(n, mu, lambda)
  } else {
    out <- rinvGAUS(n, mu, lambda)
  }
  return(out)
}


#' Truncated Normal Distribution
#'
#' Random number generation, probability density and cumulative density
#' functions for truncated normal distribution.
#'
#' @param x,q vector of quantiles;
#' @param n number of observations. n must be a scalar.
#' @param mean mean (must be scalar).
#' @param sd standard deviation (must be scalar).
#' @param lower lower truncation value (must be scalar).
#' @param upper upper truncation value (must be scalar).
#' @param lt lower tail. If TRUE (default) probabilities are \code{P[X <= x]},
#' otherwise, \code{P[X > x]}.
#' @param lp log probability. If TRUE (default is FALSE) probabilities p are
#' given as \code{log(p)}.
#' @return a column vector.
#' @examples
#' ## rtn example
#' dat1 <- rtnorm(1e5, 0, 1, 0, Inf)
#' ## dat2 <- msm::rtnorm(n, mean, sd, lower, upper)
#' ## den2 <- density(dat2)
#' hist(dat1, breaks="fd", freq=F)
#' ## lines(den2$x, den2$y, lwd=2.5)
#' ## res <- microbenchmark(
#' ##     rtn(n, mean, sd, lower, upper),
#' ##     msm::rtnorm(n, mean, sd, lower, upper))
#' ## print(res)
#'
#' ## dtn example
#' x <- seq(-5, 5, length.out=1e3)
#' dat1 <- dtnorm(x, 0, 1, -2, 2, 0)
#' ## dat2 <- msm::dtnorm(x, mean, sd, lower, upper, 0)
#' plot(x, dat1, type="l", lwd=2)
#' ## lines(x, dat2, lwd=2, lty="dashed", col="red")
#'
#' ## res <- microbenchmark(
#' ##     dtnorm(x, mean, sd, lower, upper, 0),
#' ##     msm::dtnorm(x, mean, sd, lower, upper, 0))
#' ## print(res)
#'
#' ## ptn example
#' x <- seq(-50, 10, length.out=1e3)
#' mean <- 0
#' sd <- 1
#' lower <- 0
#' upper <- 5
#' dat1 <- ptnorm(x, 0, 1, 0, 5, lp=TRUE)
#' ## dat2 <- msm::ptnorm(x, mean, sd, lower, upper, log.p=TRUE)
#' ## all(dat1[,1] == dat2)
#'
#' plot(x, log(dat1[,1]))
#' ## lines(x, log(dat2), col="red", lwd=2)
#' ## mtext("pnorm(x, log=TRUE)", adj = 0)
#' ## mtext("log(pnorm(x))", col = "red", adj = 1)
#' @export
dtnorm <- function(x, mean = 0, sd = 1, lower = -Inf, upper = Inf, log = FALSE) {
  if (any(sd < 0)) stop("sd must be greater than 0.\n")
  if (any(upper < lower)) stop("upper must be greater than lower.\n");
  if ( any(sd == -Inf) || any(sd == Inf) ) stop("sd must have a finite value.\n");
  if ( any(mean== -Inf) || any(mean == Inf) ) stop("mean must have a finite value.\n");
  
  if ( length(mean) != length(sd) || length(mean) != length(lower) || 
      length(mean) != length(upper)) {
    stop ("mean, sd, lower and upper must have same number of elements")
  }

  if (length(mean) == 1) {
    out <- dtn_(x, mean, sd, lower, upper, log)[,1]
  } else {
    out <- dtn(x, mean, sd, lower, upper, log)[,1]
  }

  return(out)
}

#' @rdname dtnorm
#' @export
ptnorm <- function(x, mean = 0, sd = 1, lower = -Inf, upper = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  if (any(sd < 0)) stop("sd must be greater than 0.\n")
  if (any(upper < lower)) stop("upper must be greater than lower.\n");
  if ( any(sd == -Inf) || any(sd == Inf) ) stop("sd must have a finite value.\n");
  if ( any(mean== -Inf) || any(mean == Inf) ) stop("mean must have a finite value.\n");
  
  if ( length(mean) != length(sd) || length(mean) != length(lower) || 
      length(mean) != length(upper)) {
    stop ("mean, sd, lower and upper must have same number of elements")
  }
  
  if (length(mean) == 1) {
    out <- ptn_(x, mean, sd, lower, upper, lower.tail, log.p)[,1]
  } else {
    out <- ptn(x, mean, sd, lower, upper, lower.tail, log.p)[,1]
  }
  
  return(out)
}
