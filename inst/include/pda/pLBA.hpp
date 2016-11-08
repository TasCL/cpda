#include <RcppArmadillo.h>

double SwitchModelPrior(Rcpp::NumericVector pVec) ;
Rcpp::List initializeStructures(const int nmc, const int npar, const int nchain) ;
arma::vec logLik_fft(Rcpp::List data, Rcpp::NumericVector pVec, Rcpp::List MCMC_params) ;
arma::mat init(Rcpp::List data, Rcpp::List MCMC_params,
  Rcpp::List Model_specifics, double gammaMult, int report) ;
