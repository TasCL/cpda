#include <cpda.hpp>

// Protect against compilers without OpenMP; e.g., OS X  clang
#ifdef _OPENMP
#include <omp.h>
#endif

//' Calculate edges for simulated histogram 
//'
//' This is an internal function. The user may not use it.
//' 
//' @param z a grid a scalar, usually created by 
//' \code{z = arma::linspace<arma::vec>(z0, z1, 1<<(int)p);}
//' 
//' @export
// [[Rcpp::export]]
arma::vec getEdges(arma::vec z) {
  double halfdt = (z[1] - z[0])/2;   
  arma::vec dtVec(z.n_elem);
  dtVec.fill(halfdt);
  arma::vec term1 = z - dtVec;
  arma::vec term2(1);
  double last     = z[z.n_elem - 1] + halfdt;
  term2.fill(last);
  return arma::join_cols(term1, term2) ;
}

//' @export
// [[Rcpp::export]]
arma::vec getFilter(double m, double M, double h, double p) {
    // See gaussian filter equation in https://en.wikipedia.org/wiki/Gaussian_filter
  double N_grid  = std::pow(2.0, p);
  double tmp0    = arma::datum::pi * N_grid/(M-m) ;
  arma::vec tmp1 = arma::linspace<arma::vec>(0, 1, 1 + N_grid/2) ;

  arma::vec tmp0Vec(tmp1.n_elem);
  tmp0Vec.fill(tmp0);
  arma::vec freq = tmp0Vec % tmp1 ;
  arma::vec freq2= arma::pow(freq, 2.0) ; // s^2 on p17

  double h2      = std::pow(h, 2.0) ;
  arma::vec fil0 = arma::exp(-0.5 * h2 * freq2) ;
  arma::vec fil1 = arma::flipud(fil0.rows(1, (fil0.size() - 2)));

  arma::vec out  = arma::join_cols(fil0, fil1) ;
  return out ;
}

//' Compute Summed Log-Likelihood Using KDE-FFT
//'
//' This function uses KDE-FFT method to approximate probability density.
//' \code{logLik_fft} returns only one scalar value, which is summed, logged
//' likelihood.
//'
//' @param y a vector storing empirical data.
//' @param yhat a vector storing simulated data (e.g., \code{rnorm(100)}).
//' @param h the bandwidth of kernel density estimation. If not given,
//' \code{logLik_fft} will detect the default, \code{h==0} and accordingly
//' convert it to Sliverman's rule of thumb; otherwise the function uses
//' bandwidth entered by the user.
//' @param m a multiplier to adjust \code{h} proportationally. Default is 0.8.
//' If one wish not adjust bandwidth, s/he has to enter \code{m=1}.
//' @param p a precision parameter defines the number of grid as power of 2.
//' Default value is 10 (i.e., ngrid==2^10).
//' @param n the number of simulation. When n=0, where the function will count
//' the numbers of observation in the simulated histogram. If simulating 
//' a defective distribution, one should enter the total number of simulation. 
//' @return Summed, logged likelihood
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
                  double m=0.8, double p=10, int n=0) {
  if (h==0) { h = bwNRD0(yhat, m); } else { h = m*h; }
  double z0, z1, out;
  arma::vec y_, yhat_, z, filter, signal, PDF0, PDF;

  y_     = arma::sort(y) ;
  yhat_  = arma::sort(yhat) ;
  z0     = y.min() - 3 * h;
  z1     = y.max() + 3 * h;
  z      = arma::linspace<arma::vec>(z0, z1, 1<<(int)p) ;
  filter = getFilter(z0, z1, h, p) ;  // Gaussian filter in freq domain
  signal = histd(yhat_, z, n) ;
  PDF0   = arma::real(arma::ifft(filter % arma::fft(signal))) ; // smoothed
  arma::interp1(z, PDF0, y_, PDF);
  out = arma::accu(arma::log(pmax(PDF, std::pow(10, -5)))) ;
  return out;
}

//' Compute Likelihood, using FFT method
//'
//' \code{lik_fft2} uses an identical algorithm as \code{logLik_fft}, but
//' return a matrix. Differing from \code{lik_pw}, \code{lik_fft2}
//' and \code{logLik_fft},
//' \enumerate{
//' \item takes Monte Carlo simulations,
//' \item transforms them to spectral domain using FFT,
//' \item applies a standard Gaussian kernel to smooth them,
//' \item transforms them back to signal domain and
//' \item interpolates linearly the simulation histogram to the observation 
//' to obtain estimated likelihoods.
//' }
//'
//' @param y a vector storing empirical data.
//' @param yhat a vector storing simulated data (e.g., \code{rnorm(100)}).
//' @param h the bandwidth for kernel density estimation. If not given,
//' \code{Lik_fft2} will use Sliverman's rule of thumb
//' @param m a multiplier to adjust \code{h} proportationally. Default is 0.8.
//' This applies also when the user enters his/her own bandwidth.
//' @param p a precision parameter defines the number of grid as power of 2.
//' Default value is 10 (i.e., 2^10).
//' @param n the number of simulation. When n=0, where the function will count
//' the numbers of observation in the simulated histogram. If simulating 
//' a defective distribution, one should enter the simulation numbers. 
//' @return A n x 2 matrix. Keep original order without sorting.
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
//' samp <- rnorm(1e5)              ## Monte Carlo simulation
//' h <- 0.8*bw.nrd0(samp)          ## Define bandwidth using R's bw.nrd0
//'
//' ## point-wise and fft
//' pw1  <- lik_pw(x, samp, h=h, m=1)
//' fft1 <- lik_fft(x, samp, h=h, m=1)
//'
//' ## the tails may be a bit off.
//' plot(x, pw1, type="l", lty="dotted")
//' lines(x, fft1[,2], col="darkgreen", lty="dashed")
//' lines(x, dnorm(x), col="red", lwd=2)
//'
//' ###################
//' ## Example 2     ##
//' ###################
//' rm(list=ls())
//' ## Assuming that this is empirical data
//' y <- rtdists::rLBA(1e3, A=.5, b=1, t0=.25, mean_v=c(2.4, 1.2), sd_v=c(1,1))
//' rt1  <- y[y$response==1,"rt"]
//' rt2  <- y[y$response==2,"rt"]
//' srt1 <- sort(rt1)
//' srt2 <- sort(rt2)
//' 
//' xlabel <- "RT"; ylabel <- "Density"
//' par(mfrow=c(1,2))
//' hist(rt1, "fd", freq=F, xlim=c(0.1,1.5), main="Choice 1", 
//'     xlab=xlabel, ylab=ylabel)
//' hist(rt2, "fd", freq=F, xlim=c(0.1,1.5), main="Choice 2",
//'     xlab=xlabel, ylab=ylabel)
//' par(mfrow=c(1,1))
//' 
//' ## Now draw simulations from another rLBA
//' n <- 1e5
//' samp <- rtdists::rLBA(n, A=.5, b=1, t0=.25, mean_v=c(2.4, 1.2), sd_v=c(1,1))
//' samp1 <- samp[samp[,2]==1, 1]
//' samp2 <- samp[samp[,2]==2, 1]
//' fft1 <- lik_fft(srt1, samp1, n=n)[,2]
//' fft2 <- lik_fft(srt2, samp2, n=n)[,2]
//' 
//' ## Just another way to simulate. This is to compare with gpda
//' pvec <- c(b=1, A=.5, v1=2.4, v2=1.2, sv1=1, sv2=1, t0=.25)
//' samp <- cpda::rlba1(n, pvec)
//' samp1 <- samp[samp[,2]==1, 1]
//' samp2 <- samp[samp[,2]==2, 1]
//' fft3 <- lik_fft(srt1, samp1, n=n)[,2]
//' fft4 <- lik_fft(srt2, samp2, n=n)[,2]
//' 
//' ## Calculate theoretical densities
//' den0 <- rtdists::dLBA(y$rt, y$response, A=.5, b=1, t0=.25, 
//'    mean_v=c(2.4, 1.2), sd_v=c(1, 1))
//'    
//' df0 <- cbind(y, den0)
//' df1 <- df0[df0[,2]==1,]
//' df2 <- df0[df0[,2]==2,]
//' den1 <- df1[order(df1[,1]),3]
//' den2 <- df2[order(df2[,1]),3]
//' plot(srt1,  den1, type="l")
//' lines(srt2, den2)
//' lines(srt1, fft1, col="red")
//' lines(srt2, fft2, col="red")
//' lines(srt1, fft3, col="blue")
//' lines(srt2, fft4, col="blue")
//'
// [[Rcpp::export]]
arma::vec lik_fft(arma::vec y, arma::vec yhat, double h=0,
                      double m=0.8, double p=10, int n=0) {
  if (h==0) { h = bwNRD0(yhat, m); } else { h = m*h; }

  double z0, z1;
  arma::vec y_, yhat_, z, filter, signal, PDF0, PDF, PDF_, PDF_unsorted;
  arma::mat PDFMat;
  arma::uvec y_idx;
  
  // y_     = arma::sort(y) ;
  // y_idx  = arma::sort_index(y) ;
  // yhat_  = arma::sort(yhat) ;
  z0     = y.min() - 3 * h;
  z1     = y.max() + 3 * h;
  z      = arma::linspace<arma::vec>(z0, z1, 1<<(int)p) ;
  filter = getFilter(z0, z1, h, p) ;  // Gauss filter
  // signal = histd(y_, yhat_, z, n) ;
  signal = histd(yhat, z, n) ;
  PDF0   = arma::real(arma::ifft(filter % arma::fft(signal))) ; // smoothed
  
  // arma::interp1(z, PDF0, y_, PDF);
  arma::interp1(z, PDF0, y, PDF);
  PDF_ = pmax(PDF, std::pow(10, -5)) ;
  // PDF_unsorted = PDF_.elem(y_idx); // re-arrange back to unsorted sequence
  // PDFMat = arma::join_horiz(y, PDF_unsorted); // return data matching PDF
  // PDFMat = arma::join_horiz(y, PDF_); // return data matching PDF
  return PDF_;
}

inline double gaussian(double y, arma::vec yhat, double h) {
  int ns = yhat.n_elem;
  arma::vec result(ns);

  for(arma::vec::iterator it=yhat.begin(); it!=yhat.end(); ++it)
  {
    int i = std::distance(yhat.begin(), it);     // (1/h) * K(z/h); K_h(z)
    result[i] = ( (1/(sqrt(2.0*arma::datum::pi))) * exp( -pow((y - *it)/h, 2.0) / 2.0 ) ) / h;
  }
  // (1/N_s) * sigma K_h (x-x.tidle_j)
  return ( arma::accu(result) / ns);
}
inline double gaussian_omp(double y, arma::vec yhat, double h) {
  int ns = yhat.n_elem;   // standard gaussian kernel mean=0; sigma==1
  arma::vec result(ns);
  
#ifdef _OPENMP
#pragma omp parallel for default(shared) firstprivate(y, yhat, h)
#endif
  for(int j=0; j < ns; j++)
  {
    result[j] = ( (1/(sqrt(2*arma::datum::pi))) *
      exp( -pow((y - yhat[j])/h, 2) / 2 ) ) / h;
  }
  
  return ( arma::accu(result) / ns);
}

//' Point-wise Probability Density Approximation
//'
//' \code{lik_pw} takes each observation and sequentially (or concurrently)
//' via Open MP) conduct KDEs. This function is for testing purpose. If you 
//' wish to conduct KDE smoothing point-by-point, use 'density'
//' in \code{stats}, which uses a similar FFT algorithm as \code{lik_fft}. 
//'
//' @param y a vector storing empirical observations (e.g., RTs).
//' @param yhat a vector storing simulations.
//' @param h kernel bandwidth. Default value is 0.8 times Silverman's Rule of
//' Thumb, based on simulated data (i.e., yhat).
//' @param m a bandwidth multiplier. Default is 0.8.
//' @param n the number of simulation. When n=0, where the function will count
//' the numbers of observation in the simulated histogram. If simulating 
//' a defective distribution, one should enter the simulation numbers. 
//' @param parallel a switch for parallel processing via OMP. Default is FALSE. 
//' @return a vector storing log-likelihoods
//' @references 
//' Holmes, W. (2015). A practical guide to the Probability Density
//' Approximation (PDA) with improved implementation and error characterization.
//' \emph{Journal of Mathematical Psychology}, vol. 68-69, 13--24,
//' doi: \url{http://dx.doi.org/10.1016/j.jmp.2015.08.006}. \cr
//' 
//' Turner, B. M. & Sederberg, P. B. (2014). A generalized, likelihood-free 
//' method for posterior estimation. \emph{Psychonomic Bulletin Review}, 21(2), 
//' 227-250. doi: \url{http://dx.doi.org/10.3758/s13423-013-0530-0}. \cr
//' 
//' Brown, S. & Heathcote, A. (2008). The simplest complete model of choice 
//' response time: Linear ballistic accumulation. \emph{Cognitive Psychology}, 
//' 57, 153-178. doi: \url{http://dx.doi.org/10.1016/j.cogpsych.2007.12.002}.
//' @export
//' @examples
//' ## Example 1 tests if Lik_pw match 'dnorm' and shows how to adjust 
//' ## bandwidth  
//' rm(list=ls())
//' lik_pw(0, rnorm(1e6)) ## 0.3974386
//' dnorm(0)              ## 0.3989423
//'
//' h <- 0.8*bw.nrd0(rnorm(1e6));            ## 0.04542646 
//' lik_pw(0, rnorm(1e6), h=h)      ## 0.3996198
//' lik_pw(0, rnorm(1e6), h=h, m=1) ## 0.3985692
//'
//' ## Example 2 demostrates how to use Lik_pw to get pLBA likelihoods
//' data(lba)
//' str(plba)
//' ## List of 4
//' ## $ DT1 : num [1:695] 0.488 0.801 0.376 0.507 0.532 ...
//' ## $ DT2 : num [1:305] 0.538 0.77 0.568 0.271 0.881 ...
//' ## $ eDT1: num [1:7020] 0.475 0.346 0.42 0.401 0.368 ...
//' ## $ eDT2: num [1:2980] 0.703 0.693 0.704 0.462 0.468 ...
//'
//' ## Use pointwise pda to get likelihoods for each data point
//' ## This algorithm calculates via a standard gaussian kernel directly
//' ## (1) First argument, plba$DT1, is the data.
//' ## (2) Second argument, plba$eDT1, is the simulation.
//' ## (3) The output is likelihood.
//'
//' output <- lik_pw(plba$DT1, plba$eDT1)
//' sum(output) ## Get summed, logged likelihood
//'
//' #########################
//' ## Example 3           ##
//' #########################
//' rm(list=ls())
//' n      <- 1e5
//' x      <- seq(-3,3, length.out=100) ## Support
//' xlabel <- "Observations" 
//' ylabel <- "Density" 
//' 
//' ## Approximate Gaussian densities
//' samp <- rnorm(n)
//' pw1  <- lik_pw(x, samp)
//' pw2  <- approx(density(samp)$x, density(samp)$y, x)$y
//' plot(x,  pw1, type="l", lty="longdash", xlab=xlabel, ylab=ylabel)
//' lines(x, pw2, lwd=1.5, lty="dotted")
//' lines(x, dnorm(x), lwd=2)
//' 
//' samp <- gamlss.dist::rexGAUS(n, mu=-2, sigma=1, nu=1)
//' system.time(pw1 <- lik_pw(x, samp, parallel=TRUE))
//' system.time(pw2 <- approx(density(samp)$x, density(samp)$y, x)$y)
//' plot(x,  pw1, type="l", lty="longdash", xlab=xlabel, ylab=ylabel)
//' lines(x, pw2, lwd=1.5, lty="dotted")
//' lines(x, gamlss.dist::dexGAUS(x, mu=-2, sigma=1, nu=1), col="red")
//'
//' ## Approximate densities of linear regression with Gaussian noise: 
//' ## y = ax + b + N(0, s)
//' theta <- c(a=7.5, b=3.5, s=5)
//' y     <- rnorm(length(x), mean=theta[2]+theta[1]*x, sd=theta[3])
//' dat   <- cbind(x, y)
//' plot(x, y, main="Linear Regression")
//' 
//' ## -- Because means for each data point differ, we need to generate 
//' ## simulation (ie samp) for each data point. That is, the likelihood for  
//' ## each data point is calculated by different simulated likelihood functions
//' ## -- We use sapply to gain a little speedy up. 
//' ## -- 'density' function cannot calculate density for only 1 data point 
//' pw1 <- sapply(x, function(s) {
//'   samp  <- rnorm(n, mean=theta[2]+theta[1]*s, sd=theta[3])
//'   lik_pw(s, samp)
//' })
//' plot(x,  pw1,  type="l", lty="longdash", xlab=xlabel, ylab=ylabel)
//' lines(x,  dnorm(x, theta[2]+theta[1]*x, theta[3]), col="red")
//' 
//' #########################
//' ## Example 4: LBA      ##
//' #########################
//' rm(list=ls())
//' ## Assuming that this is empirical data
//' y    <- rtdists::rLBA(1e3, A=.5, b=1, t0=.25, mean_v=c(2.4, 1.2), sd_v=c(1,1.2))
//' rt1  <- y[y$response==1,"rt"]
//' rt2  <- y[y$response==2,"rt"]
//' srt1 <- sort(rt1)
//' srt2 <- sort(rt2)
//' summary(rt1); summary(rt2)
//'     
//' n <- 1e5
//' pvec <- c(b=1, A=.5, v1=2.4, v2=1.2, sv1=1, sv2=1.2, t0=.25)
//' samp <- cpda::rlba1(n, pvec)
//' samp1 <- samp[samp[,2]==1, 1]
//' samp2 <- samp[samp[,2]==2, 1]
//' pw1 <- lik_pw(srt1, samp1, n=n)
//' pw2 <- lik_pw(srt2, samp2, n=n)
//'     
//' den0 <- rtdists::dLBA(y$rt, y$response, A=.5, b=1, t0=.25, 
//'   mean_v=c(2.4, 1.2), sd_v=c(1, 1.2))
//' df0 <- cbind(y, den0)
//' df1 <- df0[df0[,2]==1,]
//' df2 <- df0[df0[,2]==2,]
//' den1 <- df1[order(df1[,1]),3]
//' den2 <- df2[order(df2[,1]),3]
//'   
//' plot(srt1,  den1, type="l")
//' lines(srt2, den2)
//' lines(srt1, pw1, col="red")
//' lines(srt2, pw2, col="red")
//'     
//' #########################
//' ## Example 5           ##
//' #########################
//' n   <- 1e5
//' x   <- seq(-3,3, length.out=100) ## Support
//' sam <- rnorm(n)
//' 
//' ## Tested on 12-core CPU
//' rbenchmark::benchmark(replications=rep(10, 3),
//'    pw       <- cpda::lik_pw(x, sam),
//'    pw_omp   <- cpda::lik_pw(x, sam, parallel = T),
//'    columns=c('test', 'elapsed', 'replications'))
//' ##                                           test elapsed replications
//' ## 1                   pw <- cpda::lik_pw(x, sam)   2.570           10
//' ## 3                   pw <- cpda::lik_pw(x, sam)   2.485           10
//' ## 5                   pw <- cpda::lik_pw(x, sam)   2.484           10
//' ## 2 pw_omp <- cpda::lik_pw(x, sam, parallel = T)   0.993           10
//' ## 4 pw_omp <- cpda::lik_pw(x, sam, parallel = T)   1.119           10
//' ## 6 pw_omp <- cpda::lik_pw(x, sam, parallel = T)   1.024           10
//' 
//' @export
// [[Rcpp::export]]
arma::vec lik_pw(arma::vec y, arma::vec yhat, double h=0, double m=0.8, 
  int n=0, bool parallel=0) 
{
  double prop;
  int nobs = y.n_elem;
  int nsim = yhat.n_elem;
  arma::vec out(nobs);
  if (h == 0) { h = bwNRD0(yhat, m); } else { h = m * h; }
  if (n > 0) { prop = (double)nsim/(double)n; }

  for(int i=0; i< nobs; i++)
  {
    if (parallel) {
      out[i] = (n == 0) ? gaussian_omp(y[i], yhat, h) : 
      prop * gaussian_omp(y[i], yhat, h);
    } else {
      out[i] = (n == 0) ? gaussian(y[i], yhat, h) : 
      prop * gaussian(y[i], yhat, h);
    }
  }
  return out;
}


// arma::mat logLik_fft2_old(arma::vec y, arma::vec yhat, double h=0,
//   double m=0.8, double p=10, int n=0) {
//   if (h==0) { h = bwNRD0(yhat, m); } else { h = m*h; }
//   double z0, z1;
//   arma::vec y_, yhat_, z, filter, signal, PDF0, PDF, logPDF;
//   arma::mat PDFMat;
//   
//   y_     = arma::sort(y) ;
//   yhat_  = arma::sort(yhat) ;
//   z0     = y_.min() - 3 * h;
//   z1     = y_.max() + 3 * h;
//   z      = arma::linspace<arma::vec>(z0, z1, 1<<(int)p) ;
//   filter = getFilter(z0, z1, h, p) ;  // Gauss filter
//   signal = histd(y_, yhat_, z, n) ;
//   PDF0   = arma::real(arma::ifft(filter % arma::fft(signal))) ; // smoothed
//   
//   arma::interp1(z, PDF0, y_, PDF);
//   logPDF = arma::log(pmax(PDF, std::pow(10, -5))) ;
//   PDFMat = arma::join_horiz(y_, logPDF); // return data matching PDF
//   
//   return PDFMat;
// }

// arma::mat logLik_fft2_old2(arma::vec y, arma::vec yhat, double h=0,
//   double m=0.8, double p=10, int n=0) {
//   if (h==0) { h = bwNRD0(yhat, m); } else { h = m*h; }
//   
//   double z0, z1;
//   arma::vec y_, yhat_, z, filter, signal, PDF0, PDF, logPDF, logPDF_unsorted;
//   arma::mat PDFMat;
//   arma::uvec y_idx;
//   
//   y_     = arma::sort(y) ;
//   y_idx  = arma::sort_index(y) ;
//   yhat_  = arma::sort(yhat) ;
//   z0     = y_.min() - 3 * h;
//   z1     = y_.max() + 3 * h;
//   z      = arma::linspace<arma::vec>(z0, z1, 1<<(int)p) ;
//   filter = getFilter(z0, z1, h, p) ;  // Gauss filter
//   signal = histd(y_, yhat_, z, n) ;
//   PDF0   = arma::real(arma::ifft(filter % arma::fft(signal))) ; // smoothed
//   
//   arma::interp1(z, PDF0, y_, PDF);
//   logPDF = arma::log(pmax(PDF, std::pow(10, -5))) ;
//   logPDF_unsorted = logPDF.elem(y_idx);
//   PDFMat = arma::join_horiz(y, logPDF_unsorted); // return data matching PDF
//   
//   return PDFMat;
// }

// extern "C" void logLik_pw(double *y_, double *yhat_, int *ny, int *ns,
//                           double *h_, double *m_, double *n_, double *out);
// 
// void logLik_pw(double *y_, double *yhat_, int *ny, int *ns,
//                double *h_, double *m_, double *n_,  double *out) {
//   arma::vec yhat = getVec(yhat_ ,ns);
//   double density, h;
//   if (*h_==0) { h = bwNRD0(yhat, *m_); } else { h = (*m_)*(*h_); }
// 
//   for(int i=0; i< *ny; i++)
//   {
//     density = (*n_ == 0) ? gaussian(y_[i], yhat, h) : 
//     (*ns / *n_) * gaussian(y_[i], yhat, h);
//     // if (*n_==0) {
//     //   density = gaussian(y_[i], yhat, h);
//     // } else {
//     //   density = (*ns / *n_) * gaussian(y_[i], yhat, h); // defective
//     // }
//     out[i] = std::log(density);
//   }
// }

// extern "C" void logLik_pw_omp(double *y_, double *yhat_, int *ny, int *ns,
//                          double *h_, double *m_, double *n_, double *out);
// 
// void logLik_pw_omp(double *y_, double *yhat_, int *ny, int *ns, double *h_,
//                    double *m_, double *n_, double *out) {
//   arma::vec yhat = getVec(yhat_, ns);
//   double density, h;
//   if (*h_==0) { h = bwNRD0(yhat, *m_); } else { h = (*m_)*(*h_); }
// 
//   for(int i=0; i< *ny; i++)
//   {
//     if (*n_==0) {
//       density = gaussian_omp(y_[i], yhat, h);
//     } else {
//       density = (*ns / *n_) * gaussian_omp(y_[i], yhat, h); // defective
//     }
//     out[i]  = std::log(density);
//   }
// }


// Rcpp::List logLik_fft_test(arma::vec y, arma::vec yhat, double h=0,
//   double m=0.8, double p=10, int n=0) {
//   if (h==0) { h = bwNRD0(yhat, m); } else { h = m*h; }
//   double z0, z1, out;
//   arma::vec y_, yhat_, z, filter, signal, PDF0, PDF;
//   
//   y_     = arma::sort(y);
//   yhat_  = arma::sort(yhat);
//   z0     = std::min(y_.min(), yhat_.min()) - 3 * h;
//   z1     = std::max(y_.max(), yhat_.max()) + 3 * h;
//   z      = arma::linspace<arma::vec>(z0, z1, 1<<(int)p);
//   filter = getFilter(z0, z1, h, p);  // Gaussian filter in freq domain
//   signal = histd(y_, yhat_, z, n);
//   PDF0   = arma::real(arma::ifft(filter % arma::fft(signal))); // smoothed
//   arma::interp1(z, PDF0, y_, PDF);
//   out = arma::accu(arma::log(pmax(PDF, std::pow(10, -5))));
// 
//   Rcpp::List PDFList = Rcpp::List::create(
//     Rcpp::Named("avant") = PDF0, 
//     Rcpp::Named("apres") = PDF,
//     Rcpp::Named("LL") = out);
//   return PDFList;
// }
// 
// std::vector<double> find_test(int x) {
//   std::vector<double> V(5);
// 
//   for(std::vector<double>::iterator it=V.begin(); it<V.end(); it++) {
//     int idx = std::distance(V.begin(), it);
//     *it = idx; 
//   }
//   
//   auto it = std::find_if(std::begin(V), std::end(V), [](int i){return i > 3;});
//   int  idx = std::distance(std::begin(V), it);
//   std::cout << "idx " << idx << std::endl;
//   return V;
// }


