#include <RcppArmadillo.h>

arma::mat rlba(int n, double b, double A, arma::vec mean_v, 
  arma::vec sd_v, double t0); 
arma::vec rlba_n1(int n, double b, double A, arma::vec mean_v, arma::vec sd_v,
  double t0);