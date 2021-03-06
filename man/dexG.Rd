% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{dexG}
\alias{dexG}
\alias{pexG}
\alias{qexG}
\alias{rexG}
\title{The ex-Gaussian Distribution}
\usage{
dexG(x, mu = 5, sigma = 1, tau = 1, log = FALSE)

pexG(q, mu = 5, sigma = 1, tau = 1, lower.tail = TRUE, log.p = FALSE)

qexG(p, mu = 5, sigma = 1, tau = 1, lower.tail = TRUE, log.p = FALSE)

rexG(n, mu = 5, sigma = 1, tau = 1)
}
\arguments{
\item{x, q}{vector of quantiles.}

\item{mu}{vector of mu parameters.}

\item{sigma}{vector of sigma parameters.}

\item{tau}{vector of tau parameters.}

\item{log, log.p}{logical; if TRUE, probabilities p are given as \code{log(p)}}

\item{lower.tail}{logical; if TRUE, (default), probabilities are 
\eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.}

\item{p}{vector of probabilities.}

\item{n}{number of observations.}
}
\value{
\code{dexG} gives the density. \code{pexG} gives the distribtuion
function, \code{qexG} gives the quantile function, and \code{rexG} generates
random deviates. \eqn{sigma < 0} or \eqn{tau < 0} is an error.
}
\description{
Density, distribution function, quantile function and random generation for 
the ex-Gaussian distribution with \eqn{mean = mu + tau} and 
\eqn{standard deviation} \eqn{ = sigma^{2} + tau^{2}}. See Table 1 in 
Cousineau, Brown, & Heathcote (2004) for the equation.
}
\examples{
x <- seq(100, 600, length.out = 1024)
den <- cpda::dexG(x, mu=300, sigma=35, tau=100)

## Use Table 3 in Cousineau, Brown, & Heathcote (2004) as an example
n <- 1024
mu <- 910.56
sigma <- 44.721
tau <- 89.443
sam <- cpda::rexG(n, mu, sigma, tau)
hist(sam, freq = FALSE, breaks = "FD", main = "ex-Gaussian Distribution", 
  xlab = "RT (ms)", ylab="Density")
  
}
