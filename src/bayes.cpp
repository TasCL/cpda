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

/* Bayesian Inference */

//' Prior Log-likelihood
//'
//' This function computes the prior log-likelihood for the given parameter set,
//' using uniform distributions.
//'
//' @param pVec a vector storing pLBA model parameters. The sequence is 
//' critical. It is, A1, A2, v1, v2, t01, t02. 
//' @return a prior log prior likelihood
//' @export
//' @examples
//' logLik_prior(pVec)
//' 
// [[Rcpp::export]]
double logLik_prior(arma::vec pVec) {
  // pVec <- c(A1=1.8, A2=1.8, v1=4.5, v2=2.8, t01=.45, t02=.45)
  double pA1,pA2, pv1, pv2, pt01, pt02,out;
  pA1   = (pVec[0] <  0 || pVec[0] > 2.7) ? 0 : 0.1;
  pA2   = (pVec[1] <  0 || pVec[1] > 2.7) ? 0 : 0.1;
  pv1  = (pVec[2] < -3 || pVec[2] > 7) ? 0 : 0.1;
  pv2  = (pVec[3] < -3 || pVec[3] > 7) ? 0 : 0.1;
  pt01  = (pVec[4] <  0 || pVec[4] > 1) ? 0 : 0.1;
  pt02  = (pVec[5] <  0 || pVec[5] > 1) ? 0 : 0.1;

  out = pA1*pA2*pv1*pv2*pt01*pt02;
  return std::log(out);
}

//' @export
// [[Rcpp::export]]
double logLik_prior_test(arma::vec pVec) {
 // pVec <- c(A1=1.8, A2=1.8, b=2.7, v1=4.5, v2=2.8, t01=.45, t02=.45)
  double pA1,pA2, pb, pv1, pv2, pt01, pt02,out;
  pA1 = (pVec[0] <  0 || pVec[0] > pVec[2]) ? 0 : 0.1;
  pA2 = (pVec[1] <  0 || pVec[1] > pVec[2]) ? 0 : 0.1;
  pb  = (pVec[2] <  0 || pVec[2] > 10) ? 0 : 0.1;
  pv1 = (pVec[3] < -3 || pVec[3] > 7) ? 0 : 0.1;
  pv2 = (pVec[4] < -3 || pVec[4] > 7) ? 0 : 0.1;
  pt01  = (pVec[5] <  0 || pVec[5] > 1) ? 0 : 0.1;
  pt02  = (pVec[6] <  0 || pVec[6] > 1) ? 0 : 0.1;
  
  out = pA1*pA2*pb*pv1*pv2*pt01*pt02;
  return std::log(out);
}
