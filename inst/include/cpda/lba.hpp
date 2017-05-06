#include <RcppArmadillo.h>

void set_seed(unsigned int seed);
arma::mat rlba(int n, arma::vec pVec);
arma::mat rplba1(int n, arma::vec pVec);
arma::mat rplba2(int n, arma::vec pVec);

