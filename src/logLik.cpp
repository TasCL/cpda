#include <cpda.hpp>

// Protect against compilers without OpenMP; e.g., OS X  clang
#ifdef _OPENMP
#include <omp.h>
#endif

arma::vec getEdges(arma::vec z) {
  // Get bin_edges
  double dt   = z[1] - z[0] ;
  arma::vec term1 = z - dt/2 ;
  arma::vec term2(1) ;
  double last = z[z.size()-1] + dt/2;
  term2.fill(last) ;
  arma::vec out = arma::join_cols(term1, term2) ;
  return out ;
}

/* set precision parameter here too */ 
arma::vec getFilter(double m, double M, double h, double p) {
  // cannoical Gaussian kernel
  double tmp0    = 2 * arma::datum::pi * (std::pow(2, p) / (M-m)) * 0.5 ;
  arma::vec tmp1 = arma::linspace<arma::vec>(0, 1, 1 + (std::pow(2, p)/2)) ;
  arma::vec freq = tmp0 * tmp1 ;
  arma::vec s2   = arma::pow(freq, 2) ; // s^2 on p17
  double h2      = std::pow(h, 2) ;
  arma::vec fil0 = arma::exp(-0.5 * h2 * s2) ;
  arma::vec fil1 = arma::flipud(fil0.rows(1, (fil0.size() - 2)));
  arma::vec out  = arma::join_cols(fil0, fil1) ;
  return out ;
}

arma::vec density(arma::vec y, arma::vec be, double dt) {
  // y is yhat; be is binEdges;
  arma::uvec hc       = arma::histc(y, be) ;
  arma::vec bincount  = arma::conv_to<arma::vec>::from(hc);
  int ns              = arma::sum(bincount);
  arma::vec PDF_hist  = bincount / (dt * ns);
  arma::vec out       = PDF_hist.rows(0, (PDF_hist.size() - 2)) ;
  return out ;
}

//' Compute Summed Log-Likelihood Using KDE-FFT 
//'
//' This function uses KDE-FFT method to approximate probability density. 
//' \code{logLik_fft} returns only one scalar value, which is summed, logged 
//' likelihood. Use \code{logLik_fft2} to retrieve three additional elements: 
//' \emph{PDF}, \emph{z} and \emph{PDF_hist}. See \code{\link{logLik_fft2}} for 
//' further details.  
//'
//' @param y a scalar or vector storing empirical data.
//' @param yhat a scalar or vector storing simulated data (e.g., 
//' \code{rnorm(100)}).
//' @param h the bandwidth of kernel density estimation. If not given, 
//' \code{logLik_fft} will detect the default, \code{h==0} and accordingly 
//' convert it to Sliverman's rule of thumb; otherwise the function uses 
//' bandwidth entered by the user; 
//' @param m a multiplier to adjust \code{h} proportationally. Default is 0.8. 
//' If one wish not adjust bandwidth, s/he has to enter \code{m=1}. 
//' @param p a precision parameter defines the number of grid as power of 2.
//' Default value is 10 (i.e., ngrid==2^10). 
//' @return Summed log-likelihood
//' @references Holmes, W. (2015). A practical guide to the Probability Density
//' Approximation (PDA) with improved implementation and error characterization.
//' \emph{Journal of Mathematical Psychology}, \bold{68-69}, 13--24,
//' doi: http://dx.doi.org/10.1016/j.jmp.2015.08.006.
//' @seealso
//' \code{\link{bw.nrd}}, \code{\link{logLik_fft2}}, 
//' \code{\link{bandwidth.nrd}}.
//' 
//' @export
//' @examples
//' ## See logLik_fft2 for more examples
//' ## Use piecewise LBA data as an example
//' data(lba)
//' logLik_fft(plba$DT1, plba$eDT1)
//' logLik_fft(plba$DT2, plba$eDT2)
// [[Rcpp::export]]
double logLik_fft(arma::vec y, arma::vec yhat, double h=0, 
                  double m=0.8, double p=10) {
  
  if (h==0) { h = bwNRD0(yhat, m); } else { h = m*h; } 
  arma::vec y_     = arma::sort(y) ;
  arma::vec yhat_  = arma::sort(yhat) ;
  double z0 = std::min(y_.min(), yhat_.min()) - 3 * h;
  double z1 = std::max(y_.max(), yhat_.max()) + 3 * h;
  
  arma::vec z        = arma::linspace<arma::vec>(z0, z1, std::pow(2, p)) ;
  arma::vec binEdges = getEdges(z) ;
  arma::vec filter   = getFilter(z0, z1, h, p) ;  // Gauss filter
  double dt          = z[1] - z[0] ;
  arma::vec PDF_tmp ;
  
  arma::vec signal          = density(yhat_, binEdges, dt) ;
  arma::cx_vec PDF_fft      = arma::fft(signal) ;
  arma::cx_vec PDF_fft_filt = filter % PDF_fft ;
  arma::vec PDF_smoothed    = arma::real(arma::ifft(PDF_fft_filt)) ;
  arma::interp1(z, PDF_smoothed, y_, PDF_tmp);
  arma::vec PDF = arma::log(pmax(PDF_tmp, std::pow(10, -10))) ;
  return arma::accu(PDF);
}

//' Compute Summed Log-Likelihood and Return Likelihoods for Each Data Point  
//'
//' This function is identical to \code{logLik_fft}, except returning more
//' elements. Unlike \code{logLik_pw} calculates directly Gaussian kernel
//' (see p. 14, eq. 1.2 in Holmes, 2015), \code{logLik_fft2} and 
//' \code{logLie_fft},  
//' \enumerate{
//' \item takes Monte Carlo simulations,
//' \item transforms them to spectral domain using FFT,
//' \item applies Gaussian kernel to smooth them,
//' \item transforms them back to signal domain and
//' \item interpolates linearly the simulation density/histogram to the observed data 
//' points to obtain estimated likelihoods.
//' }
//'  
//' \code{logLik_fft2} returns four elements:
//' \itemize{
//' \item \bold{\emph{LL}}, summed, logged likelihood. This is the same as the 
//' return value from \code{logLik_fft}. 
//' \item \bold{\emph{PDF}}, a numeric vector storing logged probability 
//' densities for individual data point.
//' \item \bold{\emph{z}}, a numeric vector storing centre points of the 
//' simulated histogram (i.e., grid centre)
//' \item \bold{\emph{PDF_hist}} a numeric vector stoing the count of simulated
//' data point in each histogram bin 
//' } 
//'
//' @param y a scalar or vector storing empirical data.
//' @param yhat a scalar or vector storing simulated data (e.g., 
//' \code{rnorm(100)}).
//' @param h the bandwidth of kernel density estimation. If not given, 
//' \code{logLik_fft} will use Sliverman's rule of thumb; otherwise the 
//' function uses the input from the user; 
//' @param m a multiplier to adjust \code{h} proportationally. Default is 0.8. 
//' This applies also when one enters his/her own bandwidth. So if one wish 
//' not to adjust bandwidth, enter \code{m=1}. 
//' @param p a precision parameter defines the number of grid as power of 2.
//' Default value is 10 (i.e., 2^10). 
//' @return A list with 4 elements.
//' @references Holmes, W. (2015). A practical guide to the Probability Density
//' Approximation (PDA) with improved implementation and error characterization.
//' \emph{Journal of Mathematical Psychology}, \bold{68-69}, 13--24,
//' doi: http://dx.doi.org/10.1016/j.jmp.2015.08.006.
//' @seealso
//' \code{\link{logLik_pw}}, \code{\link{logLik_fft}},
//' \code{\link{bw.nrd}}, \code{\link{bandwidth.nrd}}.
//' 
//' @export
//' @examples
//' ###################
//' ## Example 1     ##
//' ###################
//' x <- seq(-3, 3, length.out=100) ## Data
//' samp <- rnorm(1e6)              ## Monte Carlo simulation 
//' h <- 0.8*bw.nrd0(samp)          ## Define bandwidth using R's bw.nrd0
//' 
//' ## First, I demonstrate how to use point-wise logLik
//' ## Note logLik_pw returns log-likelihood.
//' system.time(pw1  <- logLik_pw(x, samp, h=h, m=1))
//' ##   user  system elapsed 
//' ##  3.480   0.120   3.605 
//' 
//' ## Second, I demonstrate how to use KDE-FFT. logLik_fft2 retrieves 
//' ## log-likelihoods for all data point  
//' system.time(fft1 <- logLik_fft2(x, samp, h=h, m=1)[["PDF"]])
//' ##   user  system elapsed 
//' ##  1.056   0.024   1.079 
//' 
//' ## Third, I demonstrate how to build-in routine binding logLik_fft and 
//' ## rnorm to get Gaussian density. Enter ?logLik_norm2 to see help page
//' ## pVec stands for parameter vector. For a Gaussian model, pVec stores mean
//' ## and standard deviation.
//' system.time(fft2 <- logLik_norm2(x, pVec=c(0,1), 1e6, h=h, m=1)[["PDF"]])
//' ##   user  system elapsed 
//' ##  1.304   0.040   1.349  
//' 
//' ## You should get all 4 lines overlaping one another
//' ## the tails may be a bit off.
//' plot(x, pw1,type="l", lty="dotted")
//' lines(x, fft1, col="darkgreen", lty="dashed")
//' lines(x, fft2, col="blue", lty="dotdash")
//' lines(x, dnorm(x, log=TRUE), col="red")
//' 
//' ###################
//' ## Example 2     ##
//' ###################
//' ## Use piecewise LBA data as an example
//' data(lba)
//' logLik_fft(plba$DT1, plba$eDT1)
//' logLik_fft(plba$DT2, plba$eDT2)
//' plbaEg1 <- logLik_fft2(plba$DT1, plba$eDT1)
//' plbaEg2 <- logLik_fft2(plba$DT2, plba$eDT2)
//'
//' str(plbaEg1)
//' ## List of 4
//' ## $ LL      : num 280
//' ## $ PDF     : num [1:695, 1] 0.918 -0.456 0.668 0.872 0.788 ...
//' ## $ z       : num [1:1024, 1] 0.128 0.129 0.13 0.131 0.133 ...
//' ## $ PDF_hist: num [1:1024, 1] 0 0 0 0 0 0 0 0 0 0 ...
//' str(plbaEg2)
//' ## List of 4
//' ## $ LL      : num 45.3
//' ## $ PDF     : num [1:305, 1] 0.6396 -0.1174 0.6178 -1.3723 -0.0273 ...
//' ## $ z       : num [1:1024, 1] 0.121 0.124 0.126 0.129 0.131 ...
//' ## $ PDF_hist: num [1:1024, 1] 0 0 0 0 0 0 0 0 0 0 ...
// [[Rcpp::export]]
Rcpp::List logLik_fft2(arma::vec y, arma::vec yhat, double h=0, 
                       double m=0.8, double p=10) {
  if (h==0) { h = bwNRD0(yhat, m); } else { h = m*h; } 
  
  arma::vec y_     = arma::sort(y) ;
  arma::vec yhat_  = arma::sort(yhat) ;

  double z0   = std::min(y_.min(), yhat_.min()) - 3 * h;
  double z1   = std::max(y_.max(), yhat_.max()) + 3 * h;
  arma::vec z = arma::linspace<arma::vec>(z0, z1, std::pow(2, p)) ;
  arma::vec binEdges = getEdges(z) ;
  arma::vec filter   = getFilter(z0, z1, h, p) ;  // Gauss filter
  double dt          = z[1] - z[0] ;
  arma::vec signal   = density(yhat_, binEdges, dt) ;
  arma::cx_vec PDF_fft      = arma::fft(signal) ;
  arma::cx_vec PDF_fft_filt = filter % PDF_fft ;
  arma::vec PDF_smoothed    = arma::real(arma::ifft(PDF_fft_filt)) ;
  
  arma::vec PDF_tmp;   // Interpolate the grid likelihood to the data
  arma::interp1(z, PDF_smoothed, y_, PDF_tmp);
  arma::vec PDF    = arma::log(pmax(PDF_tmp, std::pow(10, -10))) ;
  arma::mat PDFMat = arma::join_horiz(y_, PDF); // return data matching PDF
  double LL        = arma::accu(PDF);

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("LL")        = LL,
    Rcpp::Named("PDF")       = PDFMat,
    Rcpp::Named("z")         = z,
    Rcpp::Named("PDF_hist")  = signal);
  return out ;
}

extern "C" void logLik_pw(double *y_, double *yhat_, int *ny, int *ns, 
                          double *h_, double *out);

void logLik_pw(double *y_, double *yhat_, int *ny, int *ns, double *h_, 
               double *out) {
  arma::vec yhat = getVec(yhat_ ,ns);
  double density;
  for(int i=0; i< *ny; i++)
  { // add std::log
    density = gaussian(y_[i], yhat, *h_);
    out[i] = std::log(density);
  }
}

extern "C" void logLik_pw_omp(double *y_, double *yhat_, int *ny, int *ns, 
                         double *h_, double *out);

void logLik_pw_omp(double *y_, double *yhat_, int *ny, int *ns, double *h_, 
               double *out) {
  arma::vec yhat = getVec(yhat_ ,ns);
  double density;
  
  for(int i=0; i< *ny; i++)
  { // add std::log
    density = gaussian_omp(y_[i], yhat, *h_);
    out[i]  = std::log(density);
  }
  
}
