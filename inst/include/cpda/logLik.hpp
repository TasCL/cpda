#include <RcppArmadillo.h>

arma::vec getEdges(arma::vec z) ;
arma::vec getFilter(double m, double M, double h, double p) ;
arma::vec density(arma::vec y, arma::vec be, double dt) ;
double logLik_fft(arma::vec y, arma::vec yhat, double h, double m, double p) ;
Rcpp::List logLik_fft2(arma::vec y, arma::vec yhat, double h, double m, double p) ;

