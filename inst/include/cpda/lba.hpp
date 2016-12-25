#include <RcppArmadillo.h>

arma::vec getEdges(arma::vec z) ;
arma::vec getFilter(double m, double M, double h) ;
arma::vec density(arma::vec y, arma::vec be, double dt) ;
arma::mat rlba(int n, arma::vec pVec);
arma::mat rplba1(int n, arma::vec pVec);
arma::mat rplba2(int n, arma::vec pVec);

