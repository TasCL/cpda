#include <cpda.hpp>
#include <boost/math/distributions/inverse_gaussian.hpp>

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector timesTwo(Rcpp::NumericVector x) {
  return x * 2;
}
