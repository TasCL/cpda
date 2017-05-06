#include <RcppArmadillo.h>

arma::vec getEdges(arma::vec z) ;
arma::vec getFilter(double m, double M, double h, double p) ;
double logLik_fft(arma::vec y, arma::vec yhat, double h, double m, 
  double p, int n) ;
arma::vec lik_fft(arma::vec y, arma::vec yhat, double h, double m, 
  double p, int n) ;
arma::vec lik_pw(arma::vec y, arma::vec yhat, double h, double m, 
  int n, bool parallel) ;