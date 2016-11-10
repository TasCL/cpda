#include <RcppArmadillo.h>

inline arma::vec getEdges(arma::vec z) ;
arma::vec getFilter(double m, double M, double bandwidth) ;
arma::vec pmax(arma::vec v, double min) ;
Rcpp::List rplba(arma::mat data, arma::vec pVec, arma::vec setting) ;
double logLik_fft(arma::vec y, arma::vec yhat, double m, double M,
  double h, int ns) ;
Rcpp::List logLik_fft2(arma::vec y, arma::vec yhat, double m, double M,
  double h, int ns) ;
double logLik_norm(arma::vec object, arma::vec pVec, arma::vec setting) ;
Rcpp::List logLik_norm2(arma::vec object, arma::vec pVec, arma::vec setting) ;

