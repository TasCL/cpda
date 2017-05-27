/*   Copyright (C) <2017>  <Yi-Shin Lin>
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

/* DE-MC algorithm */

//' Generate a Gamma Vector 
//'
//' This is part of DEMC algorithm. \code{gammavec} generates a gamma vector 
//' for element-wise compuation in Armadillo C++.
//' 
//' @param npar number of parameters. 
//' @param gammamult a tuning parameter for gamma mutation. Default value is 
//' 2.38.
//' @return a column vector
//' @export
//' @examples
//' pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51,  
//'           muw2=3.69, t_delay=0.31,  sv=1, swt=0.5)
//' gamma <- gammavec(length(pVec), 2.38)
// [[Rcpp::export]]
arma::vec gammavec(int npar, double gammamult=2.38) {
  arma::vec out(npar);
  for(int i=0; i<npar; i++) // d-dimension
  { 
    out[i] = (gammamult == 0) ? (0.5 * as_scalar(arma::randu(1)) + 0.5) : 
    gammamult/std::sqrt(2.0*(double)npar) ;
  }
  return out ;
}

//' Pick Other Chains Randomly
//'
//' This is part of DEMC algorithm. The function randomly samples \code{n}   
//' chains out of \code{length(chains)} chains, excluding the kth chain.
//' chains.
//'
//' @param k an integer indicating which chain is processed currently. This 
//' must be an integer within the range of 0 to \code{nchain-1} (i.e., C index).
//' No check for errorly using R index, because this is an internal function.
//' @param n the numbers of chain to sample.
//' @param chains a numeric vector of chain index, e.g., 0:24.  
//' @return a column vector
//' @keywords pickchains
//' @export
//' @examples
//' nchain <- 24
//' chainSeq <- (1:24)-1 ## Convert to C index
//' 
//' ## Current processing chain is the fourth chain (index=3) 
//' ## We wish to pick 2 chains out of 24 chains
//' pickchains(3, 2, chainSeq) 
//' 
// [[Rcpp::export]]
arma::vec pickchains(int k, int n, std::vector<int> chains) {
  chains.erase(chains.begin()+k) ;
  arma::vec shuffledChains = arma::shuffle(arma::conv_to<arma::vec>::from(chains)) ;
  return shuffledChains.rows(0, n-1) ;
}

//' Crossover and Migration Samplers    
//'
//' This is part of DE-MC algorithm (crossover) and Distributed Evolutionary 
//' Monte Carlo (migration). Both samplers propose a set of new parameters. 
//' 
//' @param theta a npar x nchain matrix 
//' @param gamma a gamma vector. 
//' @param k current chain index. Note chain index starts from 0.
//' @param rp ter Braak's (2006) b, setting the range of the uniform 
//' distribution that derives epsilon. This is set at 1e-3.
//' @return a column vector
//' @references ter Braak, C. J. F. (2006). A Markov Chain Monte Carlo version 
//' of the genetic algorithm Differerntial Evolution: Easy Bayesian computing
//' for real parameter spaces. Statistics and Computing, 16, 239–249. \cr\cr
//' Hu, B., & Tsui, K.-W. (2005). Distributed evolutionary Monte Carlo with 
//' applications to Bayesian analysis Technical Report Number 1112.
//' @export
//' @examples
//' init <- function(nc, npar) {
//'   theta0  <- array(dim = c(nc, npar))  
//'   theta0[,2] <- runif(nc, 0, 10);   # b
//'   theta0[,1] <- runif(nc, 0, theta0[,2]); # A
//'   theta0[,3] <- runif(nc, 0, 7);   # v1
//'   theta0[,4] <- runif(nc, 0, 7);   # v2
//'   theta0[,5] <- runif(nc, 0, .5)   # t0           
//'   return(theta0)
//' }
//' 
//' npar <- 5 
//' nc   <- npar * 1 
//' theta0  <- init(nc, npar)
//' gamma <- cpda::gammavec(npar, 2.38)
//'   
//' ## Current chain index is 0 (C-based index)
//' crossover(theta0, gamma, 0)
//' crossover(theta0, gamma, 1)
//' crossover(theta0, gamma, 2)
//' crossover(theta0, gamma, 3)
//' crossover(theta0, gamma, 4)
//' 
//' migration(theta0, gamma, 4)
//' get_subchains(nc)
//' 
// [[Rcpp::export]]
arma::mat crossover(arma::mat theta, double gammamult=2.38, double rp=0.001) {
  arma::mat theta0 = arma::trans(theta);              // theta:  nc x npar
  int npar=theta0.n_rows, nchain=theta0.n_cols; // theta0: npar x nc
  arma::vec tune = gammavec(npar, gammamult);
  arma::mat out(npar, nchain); out.fill(NA_REAL);
  std::vector<int> chains = shuffle(nchain) ; // make sure x_R0 is also random
  arma::vec subchains;

  for(int k=0; k<nchain; k++)
  {
    subchains = pickchains(k, 2, chains);
    out.col(k) = theta0.col(k) + tune % (theta0.col(subchains[0]) -
      theta0.col(subchains[1])) + (2*rp * (arma::randu(npar) - rp));
  }
  return arma::trans(out);
}


//' @rdname crossover
//' @export
// [[Rcpp::export]]
arma::vec getsubchains (int nchain) {
  // Migration algroithm - two-step shuffling
  // First shuffle decides how many chains (n) to process 
  // Second shuffle decides a subset of n chain(s) 
  arma::vec chainsSeq1 = arma::conv_to<arma::colvec>::from(shuffle(nchain));
  arma::vec chainsSeq2 = arma::conv_to<arma::colvec>::from(shuffle(nchain));
  arma::vec subchains  = chainsSeq2.rows(0, chainsSeq1[0]-1);
  return arma::sort(subchains) ;
} 


//' @rdname crossover
//' @export
// [[Rcpp::export]]
arma::mat migration(arma::mat theta, double gammamult=2.38, double rp=0.001) {
  arma::mat theta0 = arma::trans(theta);              // theta:  nc x npar
  int npar=theta0.n_rows, nchain=theta0.n_cols; // theta0: npar x nc  
  arma::vec tune = gammavec(npar, gammamult);
  arma::mat out(npar, nchain); out.fill(NA_REAL);   

  // Step 1: select a number l uniformly from 1 and k to be the number of 
  // subpopulations for migration
  arma::vec subchains = getsubchains(nchain);   
  // Connect the last number to the first number to be the migration direction

  return out;
}


/*
arma::mat migrate_primitive(arma::mat& useTheta,
  arma::vec& useLogPrior, arma::vec& useLogLike,
  Rcpp::List& pPrior, Rcpp::List& data, double rp) {
  int nChains = useTheta.n_rows ;
  int npars   = useTheta.n_cols ;
  std::vector<int> subchains    = get_subchains(nChains) ; // C index
  int nSubchains                = subchains.size() ;
  
  arma::mat thetaSet(nSubchains, npars) ;
  arma::vec currentLogPrior(nSubchains) ;
  arma::vec proposeLogPrior(nSubchains) ;
  arma::vec currentLogLike(nSubchains) ;
  arma::vec proposeLogLike(nSubchains) ;
  
  arma::vec currwLogLike(nSubchains) ;
  arma::vec propwLogLike(nSubchains) ;
  arma::vec currwLogPrior(nSubchains) ;
  arma::vec propwLogPrior(nSubchains) ;
  
  Rcpp::NumericVector model       = data.attr("model") ;
  Rcpp::NumericVector pVecNA      = model.attr("p.vector") ;
  thetaSet.fill(NA_REAL) ;
  
  for(int i=0; i<nSubchains; i++)
  {
    arma::rowvec perturbation  = Rcpp::runif(npars, -rp, rp) ;
    // create a set of particles to swap
    int ii              = subchains[i] ; // ii is non-continuous chain index
    thetaSet.row(i)     = useTheta.row(ii) + perturbation ;
    currentLogPrior(i)  = useLogPrior(ii) ; // nChains x 1
    currentLogLike(i)   = useLogLike(ii) ;
    arma::vec pVec      = vectorise(thetaSet.row(i)) ; // Proposed new prior
    proposeLogPrior(i)  = summed_log_prior(pVec, pPrior) ;
    proposeLogLike(i)   = summed_log_likelihood(pVec, data) ;
    
    propwLogPrior(i) = proposeLogPrior(i) ;
    propwLogLike(i)  = proposeLogLike(i) ;
    
    currwLogPrior(i) = currentLogPrior(i) ;
    currwLogLike(i)  = currentLogLike(i) ;
  }
  
  // ppLogLike stands for "proposed posterior log likelihood".
  // cpLogLike stands for "current  posterior log likelihood". 
  double ppLogLike = proposeLogLike(nSubchains-1) + proposeLogPrior(nSubchains-1) ;
  double cpLogLike = currentLogLike(0) + currentLogPrior(0) ;
  if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
  double rho = exp(ppLogLike - cpLogLike) ;
  
  if (rho > Rf_runif(0, 1)) {
    useTheta.row(subchains[0]) = thetaSet.row(nSubchains-1) ;
    useLogPrior(subchains[0])  = proposeLogPrior(nSubchains-1) ;
    useLogLike(subchains[0])   = proposeLogLike(nSubchains-1) ;
  }
  
  if ( nSubchains != 1 ) {
    for(int k=0; k<(nSubchains-2); k++)  // i & j were used before
    {
      // If the current selected chain is more probable than its follower,
      // replace its follower with the current chain.
      double ppLogLike = proposeLogLike(k) + proposeLogPrior(k) ;
      double cpLogLike = currentLogLike(k+1) + currentLogPrior(k+1) ;
      if (std::isnan(ppLogLike)) { ppLogLike = -INFINITY ; }
      double rho = exp (ppLogLike - cpLogLike) ;
      
      if ( rho > Rf_runif(0, 1) )
      {
        useTheta.row(subchains[k+1]) = thetaSet.row(k) ;
        useLogPrior(subchains[k+1])  = proposeLogPrior(k) ;
        useLogLike(subchains[k+1])   = proposeLogLike(k) ;
      }
    }
  }
  
  arma::mat out (nChains, 2+npars);
  out.col(0) = useLogPrior ;
  out.col(1) = useLogLike ;
  out.cols(2, 1+npars) = useTheta ;
  
  return out ;
} ;

*/
