#include <RcppArmadillo.h>

arma::vec pmax(arma::vec v, double min);
arma::vec getVec(double *x, int *nx);
double cquantile(arma::vec y, double q);
double bwNRD0(arma::vec y, double m);
std::vector<int> generateIntVec (int n);
std::vector<int> shuffle(int n);

arma::uvec histc(arma::vec x, arma::vec edge);
arma::vec histd(arma::vec yhat, arma::vec z, int n);

arma::vec getEdges(arma::vec z) ;
arma::vec getFilter(double m, double M, double h, double p) ;

