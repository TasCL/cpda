#include <pda.hpp>

//' Switch Model's Prior Distributions
//'
//' This Rcpp function computes the prior for the given parameter set. The
//' prior for each parameter here is just a uniform distribution for
//' simplicity.
//'
//' @param pVec a parameter vector. Similar to DMC's p.vector
//' @return a prior probability density summed across parameters
//' @import Rcpp
//' @export
//' @examples
//' p.vector <- c(1.51, 3.32, 1.51, 2.24, 3.69, 0.31, 0.08)
//' SwitchModelPrior(pVec=p.vector)
//' require(rbenchmark)
//' within(benchmark(R=SwitchModel_Prior(p.vector),
//'                  Cpp=SwitchModelPrior(p.vector),
//'                  replications=rep(1e4, 3),
//'                  columns=c('test', 'replications', 'elapsed'),
//'                  order=c('test', 'replications')),
//'                  { average = elapsed/replications })
//' #   test replications elapsed  average
//' # 2  Cpp        10000   0.048 4.80e-06
//' # 4  Cpp        10000   0.032 3.20e-06
//' # 6  Cpp        10000   0.056 5.60e-06
//' # 1    R        10000   0.151 1.51e-05
//' # 3    R        10000   0.182 1.82e-05
//' # 5    R        10000   0.142 1.42e-05
// [[Rcpp::export]]
double SwitchModelPrior(Rcpp::NumericVector pVec) {
   // b    A       muv1    muw1    muv2    muw2    t_delay t_ND
   // 2.7  pVec[0] pVec[1] pVec[2] pVec[3] pVec[4] pVec[5] pVec[6]
   double b = 2.7 ;
   bool A_range    = (pVec[0] <  0 || pVec[0] > 10) ;
   bool b_range    = (b <  0 || b > 10) ;
   bool muv1_range = (pVec[1] < -3 || pVec[1] > 7) ;
   bool muw1_range = (pVec[2] < -3 || pVec[2] > 7) ;
   bool muv2_range = (pVec[3] < -3 || pVec[3] > 7) ;
   bool muw2_range = (pVec[4] < -3 || pVec[4] > 7) ;
   bool t_delay_range = (pVec[5] < 0 || pVec[5] > 1) ;
   bool t_ND_range    = (pVec[6] < 0 || pVec[6] > 1) ;

   double pA = A_range ? 0 : 0.1 ;
   double pb = b_range ? 0 : 0.1 ;
   double pV1 = muv1_range ? 0 : 0.1 ;
   double pV2 = muv2_range ? 0 : 0.1 ;
   double pW1 = muw1_range ? 0 : 0.1 ;
   double pW2 = muw2_range ? 0 : 0.1 ;
   double pT = t_delay_range ? 0 : 1 ;
   double pND = t_ND_range ? 0 : 1 ;

   double prior = pA*pb*pT*pV1*pV2*pW1*pW2*pND ;
   return(prior) ;
}

//' Initialize a DMC Sample
//'
//' This functions initializes the initial conditions and various data structures
//' for the MCMC, using C++
//'
//' @param nmc number of MCMC iteration (steps).
//' @param npar number of parameter.
//' @param nchain number of chains
//' @keywords initialize_structures
//' @return a list with 7 elements: param_old, param_chain, proposal, direction,
//' LL_keep, nmc, and nchain
//' @import Rcpp
//' @export
//' @examples
//' initializeStructures(nmc=20, npar=7, nchain=3)
// [[Rcpp::export]]
Rcpp::List initializeStructures(const int nmc,
                               const int npar,
                               const int nchain)
{
  double b = 2.7 ;
  Rcpp::NumericMatrix param_old(npar, nchain) ;

  param_old(0, Rcpp::_) = b*Rcpp::rnorm(nchain) ;  // A
  param_old(1, Rcpp::_) = 5*Rcpp::runif(nchain) ;  // muv1
  param_old(2, Rcpp::_) = 5*Rcpp::runif(nchain) ;  // muw1
  param_old(3, Rcpp::_) = 5*Rcpp::runif(nchain) ;  // muv2
  param_old(4, Rcpp::_) = 5*Rcpp::runif(nchain) ;  // muw2
  param_old(5, Rcpp::_) = .5*Rcpp::runif(nchain) ; // t_delay
  param_old(6, Rcpp::_) = .5*Rcpp::runif(nchain) ; // t_nd

  arma::cube param_chain(npar, nchain, nmc) ; // storing the MCMC iteration
  arma::mat proposal(npar, nchain) ;  // storing the temporary proposal.
  arma::mat direction(npar, nchain) ;
  arma::mat LL_keep(nmc, nchain) ;

  param_chain.fill(NA_REAL) ;
  proposal.fill(NA_REAL) ;
  direction.fill(NA_REAL) ;
  LL_keep.fill(NA_REAL) ;

  // Populate first value of the chain with initial guess
  param_chain.slice(0) = Rcpp::as<arma::mat>(param_old) ;

  Rcpp::List samples_out     = Rcpp::List::create(
  Rcpp::Named("param_old")   = param_old,
  Rcpp::Named("param_chain") = param_chain,
  Rcpp::Named("proposal")    = proposal,
  Rcpp::Named("direction")   = direction,
  Rcpp::Named("LL_keep")     = LL_keep,
  Rcpp::Named("nmc")         = nmc,
  Rcpp::Named("nchains")     = nchain) ;

  return samples_out;
}

//' Retrieve Empirical and Modelled Decision Times
//'
//' Get decision times for two-accumulator LBA model
//'
//' @param data a Rcpp List for subject data
//' @param pVec a Rcpp NumvericVector for parameter vector
//' @param MCMC_params MCMC parameters
//' @return A four-element list with estimated decision time for accumulator 1,
//' accumultor 2, empirical decision time for accumualtor 1 and accumulator 2.
//' @examples
//' load("data/Data1.rda")
//' p.vector <- c(1.51, 3.32, 1.51, 2.24, 3.69, 0.31, 0.08)
//' mcmcParams <- list(sigma_exact=1, bandwidth=.02, LL_NSAMPLE=1e2,
//'                     Nstep=30, Nchain=24, noise_size=.001, burnin=10,
//'                     resample_mod=3)
//'
//' DTs <- getDTs(d, p.vector, mcmcParams)
//' str(DTs)
//' @export
// [[Rcpp::export]]
Rcpp::List getDTs(Rcpp::List data, Rcpp::NumericVector pVec,
                          Rcpp::List MCMC_params) {
  double sigma     = MCMC_params["sigma_exact"] ; // the fixed parameter sigma
  int nsample      = MCMC_params["LL_NSAMPLE"] ;

  // pVec
  // b    A       muv1    muw1    muv2    muw2    t_delay t_ND
  // 2.7  pVec[0] pVec[1] pVec[2] pVec[3] pVec[4] pVec[5] pVec[6]
  double t_ND = pVec[6] ;
  double ST = data.attr("SwitchTime") ;
  arma::vec C1time = data.attr("C1time") ;
  arma::vec C2time = data.attr("C2time") ;
  arma::vec DT1 = C1time - t_ND ; // Subtract off the non decision time
  arma::vec DT2 = C2time - t_ND ;
  double T0 = ST + pVec[5] ;  // the delay time + switch time.

  // Do sampling. Run LL_NSAMPLE sets of accumulators forward in time
  // and see which one terminates first.
  arma::vec x1 = pVec[0]*Rcpp::runif(nsample) ; // Uniform, start point
  arma::vec x2 = pVec[0]*Rcpp::runif(nsample) ; // Uniform, start point

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

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("eDT1") = eDT1,  // A1
    Rcpp::Named("eDT2") = eDT2,  // A2
    Rcpp::Named("DT1")  = DT1,   // DT1 piece 1
    Rcpp::Named("DT2")  = DT2) ; // DT2 piece 2
  return out;
}

inline arma::vec getEdges(arma::vec z) {
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

  // fliplr(): generate a copy of matrix X, with the order of the columns
  // reversed
  arma::vec freq2   = arma::pow(freq, 2) ;
  double bandwidth2 = std::pow(bandwidth, 2) ;
  arma::vec fil0    = arma::exp(-0.5 * bandwidth2 * freq2) ;
  arma::vec fil1    = flipud(fil0.rows(1, (fil0.size() - 2)));
  arma::vec out     = join_cols(fil0, fil1) ;

  return out ;
}

arma::vec pmax(arma::vec v, double min) {
  for(int i=0; i<v.size(); i++) {
    if (v[i] < min) {v[i] = min;}
  }
  return v ;
}

//' Compute Log-likelihood Using Fast Fourier Transform
//'
//' This function implements Holmes's (2015) KDE-FFT method to calculate
//' approximate probability density.
//'
//' @param DT a vector of empirical decision times
//' @param eDT a vector of modelled decision times
//' @param m min value
//' @param M max value
//' @param h KDE bandwidth
//' @param ns number of binned samples
//' @return Log-likelihood
//' @references Holmes, W. (2015). A practical guide to the Probability Density
//' Approximation (PDA) with improved implementation and error characterization.
//' \emph{Journal of Mathematical Psychology}, \bold{68-69}, 13--24,
//' doi: http://dx.doi.org/10.1016/j.jmp.2015.08.006.
//' @export
//' @examples
//' DT1  <- read.csv("data/DT1.csv", header=F)
//' DT2  <- read.csv("data/DT2.csv", header=F)
//' eDT1 <- read.csv("data/eDT1.csv", header=F)
//' eDT2 <- read.csv("data/eDT2.csv", header=F)
//' bandwidth <- .02;
//' nsample <- 1e4
//' m <- min(c(DT1$V1, DT2$V1)) - 3 * bandwidth
//' M <- max(c(DT1$V1, DT2$V1)) + 3 * bandwidth
//' logLik_fft(DT1$V1, eDT1$V1, m, M, bandwidth, nsample)
//' logLik_fft(DT2$V1, eDT2$V1, m, M, bandwidth, nsample)
// [[Rcpp::export]]
double logLik_fft(arma::vec DT, arma::vec eDT, double m, double M,
  double h, int ns) {
  arma::vec z = arma::linspace<arma::vec>(m, M, std::pow(2, 10)) ;
  arma::vec bin_edges = getEdges(z) ;
  arma::vec filter    = getFilter(m, M, h) ;  // Gauss filter
  double dt           = z[1] - z[0] ;
  arma::uvec hc       = arma::histc(eDT, bin_edges) ;
  arma::vec bincount  = arma::conv_to<arma::vec>::from(hc);
  arma::vec PDF_hist  = bincount / (dt * ns);

  arma::vec tmp             = PDF_hist.rows(0, (PDF_hist.size() - 2)) ;
  arma::cx_vec PDF_fft      = arma::fft(tmp) ;
  arma::cx_vec PDF_fft_filt = filter % PDF_fft ;
  arma::vec PDF_smoothed = arma::real(arma::ifft(PDF_fft_filt)) ;
  arma::vec PDF;   // Interpolate the grid likelihood to the data
  arma::interp1(z, PDF_smoothed, DT, PDF);
  arma::vec PDF_tmp = pmax(PDF, std::pow(10, -5)) ;
  double LL = arma::accu(arma::log(PDF_tmp));
  return LL ;
}

//' Compute FFT Log-likelihood for Piece-wise LBA Model
//'
//' Use logLik_fft to get probability density for a piecewise LBA model.
//'
//' @param data subject data
//' @param pVec a parameter vector for piece-wise LBA Model
//' @param MCMC_params MCMC parameters
//' @return Log-likelihood
//' @export
//' @examples
//' mcmcParams <- list(sigma_exact=1, bandwidth=.02, LL_NSAMPLE=1e4,
//'                   Nstep=30, Nchain=24, noise_size=.001, burnin=10,
//'                   resample_mod=3)
//'
//' load("data/Data1.rda")
//' dplyr::tbl_dt(d)
//'
//' ## Source: local data table [1,000 x 3]
//' ##      Response ResponseTime Block
//' ##         (dbl)        (dbl) (dbl)
//' ##   1         0    0.5676343     1
//' ##   2         1    0.6183177     1
//' ##   3         0    0.8806298     1
//' ##   4         0    0.4563023     1
//' ##   5         1    0.8496136     1
//' ##   6         0    0.5866219     1
//' ##   7         1    0.6482302     1
//' ##   8         1    0.3510035     1
//' ##   9         0    0.6117150     1
//' ##   10        0    0.6940521     1
//' ##   ..      ...          ...   ...
//'
//' pVec <- c(1.51, 3.32, 1.51, 2.24, 3.69, 0.31, 0.08)
//' tmp0 <- Compute_log_likelihood_FFT(d, pVec, mcmcParams)
//' tmp1 <- logLik_pLBA(d, pVec, mcmcParams)
//'
//' require(rbenchmark)
//' within(benchmark(R=Compute_log_likelihood_FFT(d, pVec, mcmcParams),
//'   Cpp=logLik_pLBA(d, pVec, mcmcParams),
//'   replications=rep(1e2, 3),
//'   columns=c('test', 'replications', 'elapsed'),
//'   order=c('test', 'replications')),
//'   { average = elapsed/replications })
//' ##     test replications elapsed average
//' ##  2  Cpp          100   0.736 0.00736
//' ##  4  Cpp          100   0.737 0.00737
//' ##  6  Cpp          100   0.736 0.00736
//' ##  1    R          100   1.037 0.01037
//' ##  3    R          100   0.922 0.00922
//' ##  5    R          100   0.922 0.00922
//'
//' @export
// [[Rcpp::export]]
double logLik_pLBA(Rcpp::List data, Rcpp::NumericVector pVec,
  Rcpp::List MCMC_params)
{
  double h       = MCMC_params["bandwidth"] ;  // Extract the KDE parameters
  int ns         = MCMC_params["LL_NSAMPLE"] ;
  Rcpp::List DTs = getDTs(data, pVec, MCMC_params) ;
  arma::vec DT1  = DTs["DT1"] ;  // C1time-t0: empirical response type 1
  arma::vec DT2  = DTs["DT2"] ;  // C1time-t0: empirical response type 2
  arma::vec eDT1 = DTs["eDT1"] ; // A1 LBA model estimates
  arma::vec eDT2 = DTs["eDT2"] ; // A2

  // Extract the minimum and maximum switch times.
  double m   = std::min(DT1.min(), DT2.min()) - 3*h ;
  double M   = std::max(DT1.max(), DT2.max()) + 3*h ;
  double LL1 = logLik_fft(DT1, eDT1, m, M, h, ns) ; // Do FFT Smoothing
  double LL2 = logLik_fft(DT2, eDT2, m, M, h, ns) ;
  return LL1 + LL2 ;
}
