#include <RcppArmadillo.h>
arma::mat choiceDT(arma::mat data, arma::vec pVec);
double logLik_norm(arma::vec object, arma::vec pVec, int ns);
Rcpp::List logLik_norm2(arma::vec object, arma::vec pVec, int ns);
double logLik_plba(arma::mat object, arma::vec pVec, int ns);
double logLik_lba(arma::mat object, arma::vec pVec, int ns);
  
