#include <RcppArmadillo.h>
inline double rtnorm0(double l, double u);
inline double rtnorm1(double l, double u);
inline double rtnorm2(double l, double u);
inline double rtnorm3(double l, double u);
inline double rtn_scalar(double mean, double sd, double l, double u);
inline double dtn_scalar(double x, double mean, double sd, double lower, 
  double upper, int lp);
double ptn_scalar(double q, double mean, double sd, double lower, double upper,
  int lt, int lp);

arma::mat rlba(int n, double b, double A, arma::vec mean_v, 
  arma::vec sd_v, double t0); 
arma::vec rlba_n1(int n, double b, double A, arma::vec mean_v, arma::vec sd_v,
  double t0);