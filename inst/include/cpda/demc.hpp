#include <RcppArmadillo.h>
arma::vec gammaVec(int n, double gammaMult);
arma::vec pick2chains(int k, int n, std::vector<int> chains);
arma::vec crossover(arma::mat useTheta, arma::vec gamma, int k, double rp);
