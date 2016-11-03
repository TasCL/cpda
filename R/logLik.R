#' Calculate Gaussian Log-likelihood and Probability Density
#'
#' A \code{logLik} function to compute density and log-likelihood, using
#' probability density approximation.
#'
#' @param object a data vector to be fit. I've sorted it internally.
#' @param pVec a parameter vector for the approximated likelihood
#' @param params a list carrying kernel density estimation parameters.
#' Nsample, bandwidth, and precision (power of 2)
#' @param ... other arguments
#' @return a list containing four elements: LL = Log likelihood, PDF = PDF of
#' the data values, grid_centers  = a histogram grid, and PDF_hist = PDF on
#' that histogram grid.
#' @keywords logLik.pLBA
#' @export
#' @references Holmes, W. (2015). A practical guide to the Probability Density
#' Approximation (PDA) with improved implementation and error characterization.
#' \emph{Journal of Mathematical Psychology}, \bold{68-69}, 13--24,
#' doi: http://dx.doi.org/10.1016/j.jmp.2015.08.006.
#' @importFrom pracma fliplr interp1
#' @importFrom graphics hist
#' @importFrom stats rnorm fft
#' @examples
#' rm(list=ls())
#'
#' ## Simulate a data vector from a Gaussian distribution
#' params <- c(mu=5, sigma=1)
#' datavec <- rnorm(1000, params[1], params[2])
#'
#' ## 1. We begin by first drawing N_s = 10,000 from the underlying distribution
#' ## to generate a synthetic data set X^~
#' ns <- 10000; ## Number of synthetic data samples used for likelihood estimation.
#' iqr <- IQR(datavec, type=1);
#' S <- sd(datavec);
#'
#' mcmcParams <- list(Nsample=ns, ## Silvermans rule of thumb
#'                    bandwidth=0.9*min(iqr, S)*ns^(-0.2),
#'                    precision=10);
#'
#' class(datavec) <- "pda"
#' ll <- logLik(datavec, params, mcmcParams);
#'
#' end1 <- length(ll$PDF_hist)-1
#' ## png(file="figs/pda.png", width=800,height=600)
#' plot(ll$gridCenters, ll$PDF_hist[1:end1], type="l", lty=2,
#'      main="Normal Distribution",xlab="x",ylab="L(x|beta)")
#' lines(datavec,ll$PDF, col="red", lwd = 1.5)
#' ## dev.off()
#' str(ll)
logLik.pda <- function(object, pVec, params, ...) {

  data <- sort(object)

  sampvec <- rnorm(params[[1]], pVec[1], pVec[2])

  m <- min(data)-3*params[[2]]
  M <- max(data)+3*params[[2]]
  ngrid <- 2^params[[3]]; ## 2^n (n>8)
  gridCenters <- seq(m, M, len = ngrid);   ## Setup the histogram bins
  dt <- gridCenters[2] - gridCenters[1];

  ## bin_edges=[gridCenters-dt/2 , gridCenters(end)+dt/2];
  term1 <- gridCenters - dt/2
  term2 <- gridCenters[length(gridCenters)] + dt/2
  bin_edges <- c(term1, term2);

  if(min(term1) > min(sampvec)) {bin_edges[1] <- min(sampvec)}
  if(max(term2) < max(sampvec)) {bin_edges[length(bin_edges)] <- max(sampvec)}

  tmp <- hist(sampvec, breaks=bin_edges, plot = FALSE, right = TRUE)
  ## Pad 0s at the lower and upper boundaries
  tmp$counts[1] <- 0
  tmp$counts[length(tmp$counts)] <- 0
  bincount <- c(tmp$counts,0);

  ## Construct the Gaussian filter that will be used to filter the data.
  ## freq=2*pi*N_grid/(M-m)*0.5 * linspace(0,1,N_grid/2+1);
  term3 <- 2*pi*ngrid / (M - m) * 0.5
  term4 <- seq(0, 1, length.out=ngrid/2+1)
  freq  <- term3 * term4;  ## s^2 in page 17 before 3.1 Resampled PDA

  fil0 <- exp(-0.5 * freq^2 * params[[2]]^2)                ## filter 0
  fil1 <- fliplr( t(as.matrix(fil0[2:(length(fil0)-1)])) )  ## filter 1
  fil  <- c(fil0, fil1)                                      ## filter

  ## Transform, smooth, and transform back to compute the PDF
  ## Create a discrete histogram representation of the PDF
  PDF_hist <- (1 * bincount) / (dt * params[[1]]);  ## dt * Nsample

  PDF_fft <- fft(PDF_hist[1:(length(PDF_hist))-1])
  PDF_fft_filt <- fil * PDF_fft;
  PDF_smoothed <- fft(PDF_fft_filt, inverse=TRUE) / length(PDF_fft_filt);
  PDF <- interp1(gridCenters, as.vector(Re(PDF_smoothed)), data,
                 method="linear");
  PDF <- pmax(PDF, 10^(-4));
  ll <- sum(log(PDF));
  return( list(ll=ll , PDF=PDF , gridCenters=gridCenters, PDF_hist=PDF_hist) )
}

