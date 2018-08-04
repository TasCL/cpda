#include <RcppArmadillo.h>
inline double rtnorm0(double l, double u);
inline double rtnorm1(double l, double u);
inline double rtnorm2(double l, double u);
inline double rtnorm3(double l, double u);
double rtn_scalar(double mean,  double sd, double l, double u);
double dtn_scalar(double x, double mean, double sd, double lower,
  double upper, bool lp);
double ptn_scalar(double q, double mean, double sd, double lower, double upper,
  bool lt, bool lp);

