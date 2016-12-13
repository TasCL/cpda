#include <RcppArmadillo.h>
arma::vec pmax(arma::vec v, double min);
arma::vec getVec(double *x, int *nx);
double gaussian(double y, arma::vec yhat, double h);
double cquantile(arma::vec y, double q);
double bwNRD0(arma::vec y, double m);
std::vector<int> generateIntVec (int n);
std::vector<int> shuffle(int n);