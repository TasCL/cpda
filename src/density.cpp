#include <cpda.hpp>

//' @export
// [[Rcpp::export]]
double gaussian(double x, arma::vec xtilde, double h) {
  unsigned int N_s = xtilde.n_elem; // See p14, equation (2.1) in Holmes (2015)
  double z, res = 0;
  
  for (size_t i = 0; i < N_s; i++) {
    z = (x - xtilde(i))/h;
    res += (1/(SQRT_2PI*h)) * exp(-0.5*z*z); 
  }
  return res/N_s;
}

//' @export
// [[Rcpp::export]]
arma::vec lik_fft(arma::vec y, arma::vec yhat, double h = 0,
                  double m = 0.8, double p = 10, unsigned int n = 0) {
  double z0, z1;
  arma::vec z, filter, signal, PDF0, PDF;
  
  h      = (h == 0) ? bwNRD0(yhat, m) : h*m;
  z0     = y.min() - 3.0 * h;
  z1     = y.max() + 3.0 * h;
  z      = arma::linspace<arma::vec>(z0, z1, 1<<(int)p) ;
  filter = getFilter(z0, z1, h, p) ;  // Gauss filter
  signal = histd(yhat, z, n) ;
  PDF0   = arma::real(arma::ifft(filter % arma::fft(signal))) ; // smoothed
  arma::interp1(z, PDF0, y, PDF);
  return pmax(PDF, std::pow(10, -5)) ;
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
arma::vec lik_pw(arma::vec x, arma::vec xtilde, double h = 0, double m = 0.8, 
  unsigned int n = 0) 
{
  unsigned int N_s = xtilde.n_elem;
  arma::vec out(x.n_elem);
  h = (h == 0) ? bwNRD0(xtilde, m) : h*m;

  for (size_t i = 0; i < x.n_elem; i++) {
    double den = gaussian(x(i), xtilde, h);
    // If n == 0 (non-defective likelihood), return den; otherwise return 
    // den / den * (N_s / n)
    out(i) = (n == 0) ? den : den * ((double)N_s/(double)n);
  }
  return out;
}



//' @export
// [[Rcpp::export]]
arma::vec n1PDF(arma::vec x, int nsim, double b, double A, arma::vec mean_v, 
  arma::vec sd_v, double t0, double h_in = 0, double k = 0.09, 
  bool debug = false) {
  // rlba_n1
  arma::vec sim = rlba_n1(nsim, b, A, mean_v, sd_v, t0);
  int nx = x.n_elem;
  arma::vec out(nx);
  int nsRT0     = sim.n_elem;
  
  if (nsRT0 <= 10) {
    out = 1e-10;
  } else {
    double minRT0 = sim.min();
    double maxRT0 = sim.max();
    double h  = (h_in == 0) ? (k*arma::stddev(sim)*std::pow(nsim, -0.2)) : h_in;
    if (debug) {Rcpp::Rcout << "h: " << h << std::endl; }
    double z0 = minRT0 <= 0 ? minRT0 : minRT0 - 3.0*h; 
    if (z0 < 0) z0 = 0;
    double z1 = maxRT0 > 10.0 ? 10.0 : maxRT0 + 3.0*h;
    int ngrid = 1024;
    int half_ngrid  = 0.5*ngrid;
    arma::vec z = arma::linspace<arma::vec>(z0, z1, ngrid);
    double dt = z[1] - z[0];
    double z1minusz0 = z1 - z0;
    double fil0_constant = (-2.0*h*h*M_PI*M_PI) / (z1minusz0*z1minusz0);

    arma::vec filter0(ngrid);
    arma::vec h_binedge0(ngrid + 1);
    arma::vec signal0(ngrid);

    // Get binedge (1025) and histogram (1024) -----------------
    for(size_t i=0; i<ngrid; i++) {
      h_binedge0[i] = z0 + dt*((double)i - 0.5); // binedge
      if (i < (1 + half_ngrid)) {                // Get filter (1024)
        filter0[i] = std::exp(fil0_constant * (double)(i*i));
      } else { 
        int j = 2*(i - half_ngrid); // flipping
        filter0[i] = filter0[i-j];
      }
    }
  
    h_binedge0[ngrid] = (z0 + ((double)(ngrid - 1))*dt);
    arma::vec h_hist0 = arma::conv_to<arma::vec>::from(arma::histc(sim, h_binedge0)); // 1025
    signal0 = h_hist0.rows(0, ngrid-1) / (dt * (double)(nsim));
  
    arma::vec sPDF = arma::real(arma::ifft(filter0 % arma::fft(signal0))) ; 
    arma::vec eDen; // a container for estiamted densities
    arma::interp1(z, sPDF, x, eDen);
  
    for(size_t i=0; i< nx; i++) { 
      out[i] = (eDen[i] < 1e-10 || std::isnan(eDen[i])) ? 1e-10 : eDen[i]; 
    }
  }
  
  return out;

}

  
  
  