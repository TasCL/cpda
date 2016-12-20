/*   Copyright (C) <2016>  <Yi-Shin Lin>
 *   This program is free software; you can redistribute it and/or modify it 
 *   under the terms of the GNU General Public License as published by the Free
 *   Software Foundation; version 2
 *    
 *   This program is distributed in the hope that it will be useful, but 
 *   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
 *   for more details.
 *    
 *   You should have received a copy of the GNU General Public License along
 *   with this program; if not, write to the Free Software Foundation, Inc.,
 *   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */
#include <cpda.hpp>

/* DEMC algorithm */

//' Generate a Gamma Vector 
//'
//' This is part of DEMC algorithm. \code{gammaVec} generates a gamma vector to 
//' facilitate element-wise compuation at C++.
//' 
//' @param n number of parameters. 
//' @param gamma a tuning parameter for gamma mutation. Usually set at 2.38  
//' @return a column vector
//' @examples
//' pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51,  
//'           muw2=3.69, t_delay=0.31,  sv=1, swt=0.5)
//' gamma <- gammaVec(length(pVec), 2.38)
arma::vec gammaVec(int n, double gamma) {
  arma::vec out(n) ;     
  for(int i=0; i<n; i++) // d-dimension
  { // an uniform random (.5, 1) 
    double r = 0.5 * as_scalar(arma::randu(1)) + 0.5; 
    out[i] = std::isnan(gamma) ? r : gamma/sqrt(2*n) ;
  }
  return out ;
}

//' Pick Other Chains Randomly
//'
//' This is part of DEMC algorithm. The function randomly chooses n chains, but
//' not the current processed one (k) in a set of \code{length(chains)} 
//' chains.
//'
//' @param k an integer indicating which chain is currently running. This must
//' be an integer within the range of 0 to \code{nchain-1} (C index). No check
//' for errorly using R index, because this is an internal function.
//' @param n number of picked chains.
//' @param chains a numeric vector of chain index, e.g., 0:24.  
//' @return a column vector
//' @keywords pickchains
//' @examples
//' nchain <- 24
//' chainSeq <- (1:24)-1 ## Convert to C index
//' 
//' ## Current processing chain is the fourth chain (index=3) 
//' ## We wish to pick 2 chains out of 24 chains
//' pickchains(3, 2, chainSeq) 
//' 
arma::vec pickchains(int k, int n, std::vector<int> chains) {
  chains.erase(chains.begin()+k) ;
  arma::vec chains0  = arma::conv_to<arma::vec>::from(chains);
  arma::vec shuffledChains = arma::shuffle(chains0) ;
  arma::vec out = shuffledChains.rows(0, n-1) ;
  return out ;
}

//' Crossover Sampler   
//'
//' This is part of DEMC algorithm. \code{crossover} proposes a set of new 
//' parameters based on crossover algorithm. 
//' 
//' @param useTheta a npar x nchain matrix 
//' @param gamma a gamma vector. 
//' @param k current chain index. Note chain index starts from 0.
//' @param rp ter Braak's (2006) b, setting the range of the uniform 
//' distribution that derives epsilon. This is usually set 0.001 (Andrew 
//' Heathcote's ROT) or ter Braak suggested using 1e-4.  My experience is 1e-3
//' is OK, haven't had too many experience using 1e-4.  
//' @return a column vector
//' @references Ter Braak, C. J. F. (2006). A Markov Chain Monte Carlo version 
//' of the genetic algorithm Differerntial Evolution: easy Bayesian computing
//' for real parameter spaces. 
//' @examples
//' pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51,  
//'           muw2=3.69, t_delay=0.31, sv=1, swt=0.5)
//' setting  <- c(bandwidth=.02, ns=1e4, nmc=6, nchain=24, rp=.001, burnin=10,
//'               nthin=3, start=1, gammaMult=2.38, report=20)
//' samples  <- init(pVec, setting)
//' useTheta <- samples$theta[,,1]
//' dim(useTheta)  ## 10 x 24; npar x nchain
//' 
//' gamma <- gammaVec(length(pVec), setting[9])
//' 
//' ## Current chain index is 0 (C-based index)
//' crossover(useTheta, gamma, 0, 1e-4)
//' 
arma::vec crossover(arma::mat useTheta, arma::vec gamma, int k, double rp) {
  int npar   = useTheta.n_rows ;   // useTheta is npar x nchain; 
  int nchain = useTheta.n_cols;
  std::vector<int> chains = shuffle(nchain) ; // make sure x_R0 is also random
  arma::vec subchains = pickchains(k, 2, chains);
  // For the sake of detail balance, use k instead of subchains[0] 
  // Scheme DE1, (Storn and Price, 1995); DEMC veci/x_i (ter Braak, 2006)
  arma::vec vk = useTheta.col(k); 
  arma::vec v1 = useTheta.col(subchains[0]); // DEMC vec1/x_R1
  arma::vec v2 = useTheta.col(subchains[1]); // DEMC vec2/x_R2
  arma::vec epsilon = rp * (arma::randu(npar) - .5);
  arma::vec proposal = vk + gamma % (v1 - v2) + epsilon; // x_p
  return proposal;
}

//' Initialise a pLBA Bayeisan Sample
//'
//' This functions initialises a pLBA sample. The function is imcompleted. 
//' The user should not use.  
//'
//' @param pVec a vector storing pLBA model parameters. The sequence is 
//' critical. It is, A, b, muv1, muv2, t_ND, muw1, muw2, t_delay, sv, swt. 
//' @param setting DE and MC setting, including, nmc, nchain, etc.
//' @keywords initialize_structures
//' @return a list with 7 elements: param_old, param_chain, proposal, direction,
//' LL_keep, nmc, and nchain
//' @examples
//' pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51,  
//'           muw2=3.69, t_delay=0.31,  sv=1, swt=0.5)
//' 
//' ## Setting sequence is critical, too! 
//' setting <- c(bandwidth=.02, ns=1e5, nmc=30, nchain=24, rp=.001, burnin=10,
//' nthin=3, start=1, gammaMult=2.38, report=100)
//' tmp0 <- init(pVec, setting)
//' str(tmp0)
//' ## List of 5
//' ## $ opVec    : num [1:10, 1:24] -1.749 3.207 2.432 4.913 0.391 ...
//' ## $ theta    : num [1:10, 1:24, 1:30] -1.749 3.207 2.432 4.913 0.391 ...
//' ## $ prospoal : num [1:10, 1:24] NA NA NA NA NA NA NA NA NA NA ...
//' ## $ direction: num [1:10, 1:24] NA NA NA NA NA NA NA NA NA NA ...
//' ## $ LL_keep  : num [1:30, 1:24] NA NA NA NA NA NA NA NA NA NA ...
//' 
Rcpp::List init(arma::vec pVec, arma::vec setting)
{
  // A  b muv1  muv2  t_ND  muw1  muw2 t_delay  sv  swt
  // 0  1    2     3     4     5     6       7   8    9  
  // block_1 = c(0, 1, 2, 3, 4); ## A, b, muv1, muv2, t_ND
  // block_2 = c(5, 6, 7);       ## muw1, muw2, t_delay
  
  // h ns nmc nchain rp burnin nthin start gammaMult report
  // 0  1   2      3  4      5     6     7         8      9
  int npar   = pVec.n_elem; 
  int nchain = setting[3];
  int nmc    = setting[2]; // dmc uses nchain x npar x nmc
  
  arma::mat opVec_(nchain, npar); opVec_.fill(NA_REAL);
  opVec_.col(1)  = 5 * arma::randu(nchain);   // b
  opVec_.col(0)  = opVec_.col(1) % arma::randn(nchain) ;  // A
  opVec_.col(2)  = 5*arma::randu(nchain) ;  // muv1
  opVec_.col(5)  = 5*arma::randu(nchain) ;  // muw1
  opVec_.col(3)  = 5*arma::randu(nchain) ;  // muv2
  opVec_.col(6)  = 5*arma::randu(nchain) ;  // muw2
  opVec_.col(7)  = .5*arma::randu(nchain) ; // t_Delay
  opVec_.col(4)  = .5*arma::randu(nchain) ; // t_ND
  opVec_.col(8)  = arma::repmat(pVec.row(8), nchain, 1); // sv
  opVec_.col(9)  = arma::repmat(pVec.row(9), nchain, 1); // swt
  
  
  // Initialise the 1st theta; rprior
  arma::cube theta(npar, nchain, nmc); // param_chain (nmc x npar x nchain); 
  arma::mat theta_(npar, nchain);      // temporarily store 1st slice
  theta.fill(NA_REAL);     // Populate first value of the chain with initial guess
  theta_.row(1)  = opVec_.col(1).t(); // b
  theta_.row(0)  = opVec_.col(0).t(); // A
  theta_.row(2)  = opVec_.col(2).t(); // muv1
  theta_.row(5)  = opVec_.col(5).t(); // muw1
  theta_.row(3)  = opVec_.col(3).t(); // muv2
  theta_.row(6)  = opVec_.col(6).t(); // muw2
  theta_.row(7)  = opVec_.col(7).t(); // t_Delay
  theta_.row(4)  = opVec_.col(4).t(); // t_ND
  theta_.row(8)  = opVec_.col(8).t(); // sv
  theta_.row(9)  = opVec_.col(9).t(); // swt
  theta.slice(0) = theta_;
  arma::mat opVec = opVec_.t(); // back to npar x nchain
  
  arma::mat proposal(npar, nchain);
  arma::mat direction(npar, nchain);
  arma::mat LL_keep(nmc, nchain);      
  proposal.fill(NA_REAL);
  direction.fill(NA_REAL);
  LL_keep.fill(NA_REAL);
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("opVec")     = opVec,
    Rcpp::Named("theta")     = theta,
    Rcpp::Named("prospoal")  = proposal,
    Rcpp::Named("direction") = direction,
    Rcpp::Named("LL_keep")   = LL_keep) ;
  
  return out;
}

