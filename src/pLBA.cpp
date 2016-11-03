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

//' Compute Log-likelihood Using Fast Fourier Transform
//'
//' This function takes subject data, the parameters for the current chain
//' under consideration, and MCMC parameters and computes the log likelihood.
//'
//' @param data subject data
//' @param pVec a parameter vector for piece-wise LBA Model
//' @param MCMC_params MCMC parameters
//' @return Log-likelihood
//' @export
//' @examples
//' mcmcParams <- list(sigma_exact=1, bandwidth=.02, LL_NSAMPLE=1e5,
//'                   Nstep=30, Nchain=24, noise_size=.001, burnin=10,
//'                   resample_mod=3)
//' logLik_fft(mcmcParams)
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
//' t0 <- as.vector(attr(d, "SwitchTime"))
//' C1time <- as.vector(attr(d, "C1time"))
//' C2time <- as.vector(attr(d, "C2time"))
//' attr(d, "SwitchTime") <- NULL
//' attr(d, "C1time")     <- NULL
//' attr(d, "C2time")     <- NULL
//' str(d)
//' ## 'data.frame':	1000 obs. of  3 variables:
//' ## $ Response    : num  0 1 0 0 1 0 1 1 0 0 ...
//' ## $ ResponseTime: num  0.568 0.618 0.881 0.456 0.85 ...
//' ## $ Block       : num  1 1 1 1 1 1 1 1 1 1 ...
//' dat <- list(data=d, SwitchTime=t0, C1time=C1time, C2time=C2time)
//' p.vector <- c(1.51, 3.32, 1.51, 2.24, 3.69, 0.31, 0.08)
//' tmp0 <- logLik_fft(data=dat, pVec=p.vector, MCMC_params=mcmcParams)
//' tmp1 <- dat$C2time - p.vector[7]
//' all(tmp0==tmp2)
//' [1] TRUE
// [[Rcpp::export]]
arma::vec logLik_fft(Rcpp::List data, Rcpp::NumericVector pVec, Rcpp::List MCMC_params) {
  double sigma     = MCMC_params["sigma_exact"] ; // the fixed parameter sigma
  double bandwidth = MCMC_params["bandwidth"] ;  // Extract the KDE parameters
  int nsample      = MCMC_params["LL_NSAMPLE"] ;
  int nmc          = MCMC_params["Nstep"] ;
  int nchain       = MCMC_params["Nchain"] ;
  double rp        = MCMC_params["noise_size"] ;
  int burnin       = MCMC_params["burnin"] ;
  int resample_mod = MCMC_params["resample_mod"] ;

  // pVec
  // A muv1 muw1 muv2 muw2 t_delay t_ND
  // 0    1    2    3    4       5    6
  double t_ND = pVec[6] ;
  double T0 = data["SwitchTime"] ;
  arma::vec C1time = data["C1time"] ;
  arma::vec C2time = data["C2time"] ;
  arma::vec DT1 = C1time - t_ND ;
  arma::vec DT2 = C2time - t_ND ;

  return DT2 ;
}

//' Retrieve piecewise RTs from a pLBA Model
//'
//' Get piecewise RTs
//'
//' @param data a Rcpp List for subject data
//' @param pVec a Rcpp NumvericVector for parameter vector
//' @param MCMC_params MCMC parameters
//' @return A two-element list
//' @examples
//' load("data/Data1.rda")
//' p.vector <- c(1.51, 3.32, 1.51, 2.24, 3.69, 0.31, 0.08)
//' mcmcParams <- list(sigma_exact=1, bandwidth=.02, LL_NSAMPLE=1e2,
//'                     Nstep=30, Nchain=24, noise_size=.001, burnin=10,
//'                     resample_mod=3)
//'
//' tmp0 <- getPRT(d, p.vector, mcmcParams)
//' str(tmp0[[1]])
//' str(tmp0[[2]])
//' @export
// [[Rcpp::export]]
Rcpp::List getPRT(Rcpp::List data, Rcpp::NumericVector pVec,
                          Rcpp::List MCMC_params) {
  double sigma     = MCMC_params["sigma_exact"] ; // the fixed parameter sigma
  double bandwidth = MCMC_params["bandwidth"] ;  // Extract the KDE parameters
  int nsample      = MCMC_params["LL_NSAMPLE"] ;
  int nmc          = MCMC_params["Nstep"] ;
  int nchain       = MCMC_params["Nchain"] ;
  double rp        = MCMC_params["noise_size"] ;
  int burnin       = MCMC_params["burnin"] ;
  int thin         = MCMC_params["resample_mod"] ;

  // pVec
  // b    A       muv1    muw1    muv2    muw2    t_delay t_ND
  // 2.7  pVec[0] pVec[1] pVec[2] pVec[3] pVec[4] pVec[5] pVec[6]
  double t_ND = pVec[6] ;
  double ST = data.attr("SwitchTime") ;
  arma::vec C1time = data.attr("C1time") ;
  arma::vec C2time = data.attr("C2time") ;
  arma::vec DT1 = C1time - t_ND ; // Subtract off the non decision time
  arma::vec DT2 = C1time - t_ND ;
  double T0 = ST + pVec[5] ;  // the delay time + switch time.

  // Do sampling. Run LL_NSAMPLE sets of accumulators forward in time
  // and see which one terminates first.
  arma::vec x1_0 = pVec[0]*Rcpp::runif(nsample) ; // Uniform, start point
  arma::vec x2_0 = pVec[0]*Rcpp::runif(nsample) ; // Uniform, start point

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
  arma::vec T1_s_first  = (2.7-x1_0) / v1;
  arma::vec T2_s_first  = (2.7-x2_0) / v2;
  arma::vec T1_s_second = T0 + (2.7-x1_0-v1*T0) / w1;
  arma::vec T2_s_second = T0 + (2.7-x2_0-v2*T0) / w2;

  arma::uvec q1 = find((T1_s_first < T0) && (T1_s_first < T2_s_first)) ;
  arma::uvec q2 = find((T1_s_first > T0) && (T2_s_first > T0) &&
    (T1_s_second < T2_s_second) && (T1_s_second < 1e4)) ;

  arma::uvec q3 = find((T2_s_first < T0) && (T2_s_first < T1_s_first)) ;
  arma::uvec q4 = find((T1_s_first > T0) && (T2_s_first > T0) &&
    (T2_s_second < T1_s_second) && (T2_s_second < 1e4)) ;

  arma::vec tmp0 = T1_s_first.elem( q1 ) ;
  arma::vec tmp1 = T1_s_second.elem( q2 ) ;
  arma::vec tmp2 = T1_s_first.elem( q3 ) ;
  arma::vec tmp3 = T1_s_second.elem( q4 ) ;

  arma::vec RT_s1 = arma::join_cols(tmp0, tmp1) ;
  arma::vec RT_s2 = arma::join_cols(tmp2, tmp3) ;

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("RT_s1") = RT_s1,
    Rcpp::Named("RT_s2") = RT_s2) ;

  // Rcpp::Rcout << tmp0.size() << std::endl;
  // Rcpp::Rcout << tmp1.size() << std::endl;
  // Rcpp::Rcout << tmp2.size() << std::endl;
  // Rcpp::Rcout << tmp3.size() << std::endl;


  return out;
}

