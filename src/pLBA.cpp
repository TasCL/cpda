#include <pda.hpp>

struct genInt {
  int x ;
  genInt() {x=0;}
  int operator()() {return x++;}
} generateInteger; // class generator:

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

arma::vec getFilter(double m, double M, double bandwidth) {
  double tmp0 = 2 * arma::datum::pi * (std::pow(2, 10) / (M-m)) * 0.5 ;
  arma::vec freqSpace = arma::linspace<arma::vec>(0, 1, 1 + (std::pow(2, 10)/2)) ;
  arma::vec freq = tmp0 * freqSpace ;

  // fliplr(): generate a copy of matrix X, with reversed columns order
  arma::vec freq2   = arma::pow(freq, 2) ;
  double bandwidth2 = std::pow(bandwidth, 2) ;
  arma::vec fil0    = arma::exp(-0.5 * bandwidth2 * freq2) ;
  arma::vec fil1    = flipud(fil0.rows(1, (fil0.size() - 2)));
  arma::vec out     = join_cols(fil0, fil1) ;
  return out ;
}

arma::vec pmax(arma::vec v, double min) {
  for (arma::vec::iterator it=v.begin(); it!=v.end(); it++)
  {
    if (*it < min) *it = min ;
  }
  return v ;
}

//' Retrieve Empirical and Modelled Decision Times
//'
//' Get decision times for two-accumulator LBA model
//'
//' @param data a Rcpp List for subject data
//' @param pVec a Rcpp NumvericVector for parameter vector
//' @param setting a named numeric vector storing DE and MCMC parameters:
//' sigma_exact, bandwidth, ns (nsample), nmc, nchain, rp (DE noise,
//' usually .001), burnin, nthin, start (usually 1), gammaMult (DE gamma,
//' usually 2.38), ST (piecewise LBA switch time, depending on experimental
//' design), and report
//' @return A four-element list with estimated decision time for accumulator 1,
//' accumultor 2, empirical decision time for accumualtor 1 and accumulator 2.
//' @export
//' @examples
//' data(lba)
//' dMat <- data.matrix(d)
//'
//' ##  A       muv1    muw1    muv2    muw2    t_delay t_ND    b
//' p.vector <- c(A=1.51, muv1=3.32, muw1=1.51, muv2=2.24, muw2=3.69,
//'               t_delay=0.31, t_ND=0.08, b=2.7)
//' setting <- c(sigma_exact=1, bandwidth=.02, ns=1e2, nmc=30, nchain=24,
//'              rp=.001, burnin=10, nthin=3, start=1, gammaMult=2.38,
//'              ST=attr(d, "SwitchTime"), report=100)
//'
//' DTs <- pda::rplba(dMat, p.vector, setting)
//' str(DTs)
//'
//' ## List of 4
//' ## $ eDT1: num [1:72, 1] 0.317 0.425 0.71 0.345 0.526 ...
//' ## $ eDT2: num [1:28, 1] 0.45 0.418 0.56 0.719 0.313 ...
//' ## $ DT1 : num [1:695, 1] 0.488 0.801 0.376 0.507 0.532 ...
//' ## $ DT2 : num [1:305, 1] 0.538 0.77 0.568 0.271 0.881 ...
// [[Rcpp::export]]
Rcpp::List rplba(arma::mat data, arma::vec pVec, arma::vec setting) {
  double sigma = setting[0] ; // sigma_exact ; // the fixed parameter sigma
  int nsample  = setting[2] ; // LL_NSAMPLE ;
  arma::vec choice = data.col(0) ;   // Response ResponseTime Block; 0==error; 1==correct
  arma::vec rt     = data.col(1) ;
  arma::uvec idx0  = find(choice == 0) ; // error
  arma::uvec idx1  = find(choice == 1) ; // correct
  arma::vec C1time  = rt.rows(idx0) ; // error RT
  arma::vec C2time  = rt.rows(idx1) ; // correct RT

  // A       muv1    muw1    muv2    muw2    t_delay t_ND    b
  // pVec[0] pVec[1] pVec[2] pVec[3] pVec[4] pVec[5] pVec[6] pVec[7]
  arma::vec DT1 = C1time - pVec[6] ; // Subtract off the non decision time
  arma::vec DT2 = C2time - pVec[6] ;
  double T0 = setting[10] + pVec[5] ;  // the delay time + switch time.

  // Do sampling. Run LL_NSAMPLE sets of accumulators forward in time
  // and see which one terminates first.
  arma::vec x1 = pVec[0]*arma::randu(nsample) ; // Uniform, start point
  arma::vec x2 = pVec[0]*arma::randu(nsample) ; // Uniform, start point

  arma::vec v1(nsample);
  arma::vec w1(nsample);
  arma::vec v2(nsample);
  arma::vec w2(nsample);

  for(int i=0; i<nsample; i++) // ## Handle negative drift rates
  {
    v1[i] = rtn_scalar(pVec[1], sigma, 0, INFINITY) ;
    w1[i] = rtn_scalar(pVec[2], sigma, 0, INFINITY) ;
    v2[i] = rtn_scalar(pVec[3], sigma, 0, INFINITY) ;
    w2[i] = rtn_scalar(pVec[4], sigma, 0, INFINITY) ;
  }

  // Compute the time at which the two accumulators (v1 & v2)
  // will terminate. 2.7==b; v1 & v2 are two accumulators during 1st piece
  // ditto for w1 and w2. w1 & w2 are two accumulators during 2nd piece
  arma::vec A1P1 = (2.7-x1) / v1;
  arma::vec A2P1 = (2.7-x2) / v2;
  arma::vec A1P2 = T0 + (2.7-x1-v1*T0) / w1;
  arma::vec A2P2 = T0 + (2.7-x2-v2*T0) / w2;

  // pre switch termination
  arma::uvec A1WinB = find((A1P1 < T0) && (A1P1 < A2P1)) ; // A1 wins
  arma::uvec A2WinB = find((A2P1 < T0) && (A2P1 < A1P1)) ; // A2 wins

  // post switch termination
  arma::uvec A1WinA = find((A1P1 > T0) && (A2P1 > T0) && (A1P2 < A2P2) &&
    (A1P2 < 1e4)) ;
  arma::uvec A2WinA = find((A1P1 > T0) && (A2P1 > T0) && (A2P2 < A1P2) &&
    (A2P2 < 1e4)) ;

  arma::vec eDT1 = arma::join_cols(A1P1.elem( A1WinB ), A1P2.elem( A1WinA )) ;
  arma::vec eDT2 = arma::join_cols(A2P1.elem( A2WinB ), A2P2.elem( A2WinA )) ;

  // std::map <std::string, arma::vec> map0;
  // map0["eDT1"] = eDT1 ;
  // map0["eDT2"] = eDT2 ;
  // map0["DT1"]  = DT1 ;
  // map0["DT2"]  = DT2 ;
  // arma::vec tmp0 = map0["eDT1"] ;

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("eDT1") = eDT1,  // A1
    Rcpp::Named("eDT2") = eDT2,  // A2
    Rcpp::Named("DT1")  = DT1,   // DT1 piece 1
    Rcpp::Named("DT2")  = DT2) ; // DT2 piece 2
  return out;
}

//' Compute Log-likelihood Using KDE-based Fast Fourier Transform
//'
//' This function implements Holmes's (2015) KDE-FFT method to calculate
//' approximated probability density. The precision is set at 10 (ie 2^10).
//' Please use \code{logLik_fft2}, if you want also the PDF, grid centers, and
//' PDF_hist outputs.
//'
//' @param y a vector storing empirical data (e.g., RTs)
//' @param yhat a vector storing simulated data (e.g., simualted RTs, using a
//' LBA model).
//' @param m minimal value extracted from the empirical data vector
//' @param M max value extracted from the empirical data vector
//' @param h KDE bandwidth
//' @param ns number of simulate data
//' @return Log-likelihood
//' @references Holmes, W. (2015). A practical guide to the Probability Density
//' Approximation (PDA) with improved implementation and error characterization.
//' \emph{Journal of Mathematical Psychology}, \bold{68-69}, 13--24,
//' doi: http://dx.doi.org/10.1016/j.jmp.2015.08.006.
//' @export
//' @examples
//' ## Use piecewise LBA data as an example
//' data(lba)
//' bandwidth <- .02;
//' nsample <- 1e4
//' m <- min(c(plba$DT1, plba$DT2)) - 3 * bandwidth
//' M <- max(c(plba$DT1, plba$DT2)) + 3 * bandwidth
//' logLik_fft(plba$DT1, plba$eDT1, m, M, bandwidth, nsample)
//' logLik_fft(plba$DT2, plba$eDT2, m, M, bandwidth, nsample)
//' tmp1 <- logLik_fft2(plba$DT1, plba$eDT1, m, M, bandwidth, nsample)
//' tmp2 <- logLik_fft2(plba$DT2, plba$eDT2, m, M, bandwidth, nsample)
//' str(tmp1)
//'
//' ## List of 4
//' ## $ LL      : num 33.9
//' ## $ PDF     : num [1:695, 1] 1.761 0.445 1.368 1.681 1.542 ...
//' ## $ z       : num [1:1024, 1] 0.157 0.159 0.16 0.161 0.162 ...
//' ## $ PDF_hist: num [1:1025, 1] 0 0 0 0 0 0 0 0 0 0 ...
//' str(tmp2)
//' ## List of 4
//' ## $ LL      : num -323
//' ## $ PDF     : num [1:305, 1] 0.5526 0.2489 0.5588 0.0579 0.3016 ...
//' ## $ z       : num [1:1024, 1] 0.157 0.159 0.16 0.161 0.162 ...
//' ## $ PDF_hist: num [1:1025, 1] 0 0 0 0 0 0 0 0 0 0 ...
// [[Rcpp::export]]
double logLik_fft(arma::vec y, arma::vec yhat, double m, double M,
  double h, int ns) {
  arma::vec z = arma::linspace<arma::vec>(m, M, std::pow(2, 10)) ;
  arma::vec bin_edges = getEdges(z) ;
  arma::vec filter    = getFilter(m, M, h) ;  // Gauss filter
  double dt           = z[1] - z[0] ;
  arma::uvec hc       = arma::histc(yhat, bin_edges) ;
  arma::vec bincount  = arma::conv_to<arma::vec>::from(hc);
  arma::vec PDF_hist  = bincount / (dt * ns);

  arma::vec tmp             = PDF_hist.rows(0, (PDF_hist.size() - 2)) ;
  arma::cx_vec PDF_fft      = arma::fft(tmp) ;
  arma::cx_vec PDF_fft_filt = filter % PDF_fft ;
  arma::vec PDF_smoothed = arma::real(arma::ifft(PDF_fft_filt)) ;
  arma::vec PDF;   // Interpolate the grid likelihood to the data
  arma::interp1(z, PDF_smoothed, y, PDF);
  arma::vec PDF_tmp = pmax(PDF, std::pow(10, -5)) ;
  double LL = arma::accu(arma::log(PDF_tmp));
  return LL ;
}

//' @rdname logLik_fft
//' @export
// [[Rcpp::export]]
Rcpp::List logLik_fft2(arma::vec y, arma::vec yhat, double m, double M,
  double h, int ns) {
    arma::vec z = arma::linspace<arma::vec>(m, M, std::pow(2, 10)) ;
    arma::vec bin_edges = getEdges(z) ;
    arma::vec filter    = getFilter(m, M, h) ;  // Gauss filter
    double dt           = z[1] - z[0] ;
    arma::uvec hc       = arma::histc(yhat, bin_edges) ;
    arma::vec bincount  = arma::conv_to<arma::vec>::from(hc);
    arma::vec PDF_hist  = bincount / (dt * ns);

    arma::vec tmp             = PDF_hist.rows(0, (PDF_hist.size() - 2)) ;
    arma::cx_vec PDF_fft      = arma::fft(tmp) ;
    arma::cx_vec PDF_fft_filt = filter % PDF_fft ;
    arma::vec PDF_smoothed = arma::real(arma::ifft(PDF_fft_filt)) ;
    arma::vec PDF;   // Interpolate the grid likelihood to the data
    arma::interp1(z, PDF_smoothed, y, PDF);
    arma::vec PDF_tmp = pmax(PDF, std::pow(10, -5)) ;
    double LL = arma::accu(arma::log(PDF_tmp));

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("LL")        = LL,
        Rcpp::Named("PDF")       = PDF_tmp,
        Rcpp::Named("z")         = z,
        Rcpp::Named("PDF_hist")  = PDF_hist);
    return out ;
}

//' Compute FFT Log-likelihood for a Gaussian Distribution
//'
//' This is a wrapper function to approximate a Gaussian distribution, using
//' KDE-FFT method. To retrieve more outputs and replicate Holmes's example 1,
//' please use \code{logLik_norm2}.
//'
//' @param object a vector storing empirical data
//' @param pVec parameter vector storing mean and standard deviation
//' @param setting MCMC and DE-MC parameters
//' @return Log-likelihood; plus PDF, grid centers, and PDF_hist
//' @export
//' @examples
//' pVec <- c(mu=5, sigma=1)
//' y    <- sort(rnorm(1000, pVec[1], pVec[2]))
//' iqr  <- IQR(y, type=1)
//' S    <- sd(y)
//'
//' setting <- c(sigma_exact=1, bandwidth=.9*min(iqr,S)*1e4^(-.2), ns=1e4,
//'              nmc=30, nchain=24, rp=.001, burnin=10, nthin=3, start=1,
//'              gammaMult=2.38, ST=0, report=100)
//'
//' ll <- logLik_norm2(y, pVec, setting)
//' str(ll)
//' ## List of 4
//' ## $ LL      : num -1370
//' ## $ PDF     : num [1:1000, 1] 0.006 0.00617 0.00813 0.02605 0.03 ...
//' ## $ z       : num [1:1024, 1] 1.65 1.66 1.66 1.67 1.68 ...
//' ## $ PDF_hist: num [1:1025, 1] 0 0 0 0 0 0 0 0 0 0 ...
//'
//' end1 <- length(ll$PDF_hist)-1
//' plot(ll$z, ll$PDF_hist[1:end1], type="l", lty=2,
//' main="Normal Distribution",xlab="x",ylab="L(x|beta)")
//' lines(y, ll$PDF, col="red", lwd = 1.5)
//'
//' ##########################################################################
//'
//' logLik_norm(y, pVec, setting)
//' ## [1] -1372.969
// [[Rcpp::export]]
double logLik_norm(arma::vec object, arma::vec pVec, arma::vec setting) {
    // pVec[0] is mean, pVec[1] is sigma
    arma::vec y    = arma::sort(object) ;
    int ns         = setting[2] ;
    double h       = setting[1];
    arma::vec yhat = pVec[0]+pVec[1] * arma::randn(ns) ; // simulation
    double m  = y.min() - 3*h; // min for empirical data
    double M  = y.max() + 3*h; // max ...
    double LL = logLik_fft(y, yhat, m, M, h, ns) ;
    return LL ;
}

//' @rdname logLik_norm
//' @export
// [[Rcpp::export]]
Rcpp::List logLik_norm2(arma::vec object, arma::vec pVec, arma::vec setting) {
    arma::vec y    = arma::sort(object) ;
    int ns         = setting[2] ;
    double h       = setting[1];
    arma::vec yhat = pVec[0]+pVec[1] * arma::randn(ns) ; // simulation
    double m  = y.min() - 3*h; // min for empirical data
    double M  = y.max() + 3*h; // max ...
    Rcpp::List LL = logLik_fft2(y, yhat, m, M, h, ns) ;
    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("LL")        = LL["LL"],
      Rcpp::Named("PDF")       = LL["PDF"],
      Rcpp::Named("z")         = LL["z"], // grid centers
      Rcpp::Named("PDF_hist")  = LL["PDF_hist"]) ;
    return out ;
}
