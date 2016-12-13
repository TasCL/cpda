#include <RcppArmadillo.h>

arma::vec getEdges(arma::vec z) ;
arma::vec getFilter(double m, double M, double h) ;
arma::vec density(arma::vec y, arma::vec be, double dt) ;
arma::mat rlba(int n, arma::vec pVec);
arma::mat rplba(int n, arma::vec pVec);
arma::mat choiceDT(arma::mat data, arma::vec pVec);
double logLik_fft(arma::vec y, arma::vec yhat) ;
Rcpp::List logLik_fft2(arma::vec y, arma::vec yhat) ;

