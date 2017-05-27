#include <RcppArmadillo.h>

arma::vec lik_fft(arma::vec y, arma::vec yhat, double h, double m, double p, 
  int n);
arma::vec lik_pw(arma::vec y, arma::vec yhat, double h, double m, int n); 