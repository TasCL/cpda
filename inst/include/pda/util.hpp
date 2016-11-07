#include <RcppArmadillo.h>

inline bool CheckSimple(const double lower,  const double upper);
bool tn_check(const double mean, const double sd, const double lower,
  const double upper);
inline double alg1(const double lower, const double upper);
inline double alg2(const double lower, const double upper);
inline double alg3(const double lower, const double upper);
inline double alg4(const double lower, const double upper);
double rtn_scalar(const double mean,  const double sd, const double lower,
  const double upper);
double dtn_scalar(const double x, const double mean, const double sd,
  const double lower, const double upper, int islog);
void rtn(Rcpp::NumericVector &mean, Rcpp::NumericVector &sd,
  Rcpp::NumericVector &lower, Rcpp::NumericVector &upper,
  Rcpp::NumericVector &draws);
void dtn(Rcpp::NumericVector &x,
  Rcpp::NumericVector &mean,  Rcpp::NumericVector &sd,
  Rcpp::NumericVector &lower, Rcpp::NumericVector &upper,
  Rcpp::LogicalVector &islog, Rcpp::NumericVector &dens);
RcppExport SEXP rtn_wrapper(const SEXP mean_, const SEXP sd_,
  const SEXP lower_, const SEXP upper_);
RcppExport SEXP dtn_wrapper(const SEXP x_, const SEXP mean_, const SEXP sd_,
  const SEXP lower_, const SEXP upper_, const SEXP islog_);

Rcpp::List getPRT(Rcpp::List data, Rcpp::NumericVector pVec,
  Rcpp::List MCMC_params) ;
arma::vec getFilter(double m, double M, double bandwidth);
arma::vec logLik(Rcpp::List data, Rcpp::NumericVector pVec, Rcpp::List MCMC_params) ;
