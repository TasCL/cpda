//' Retrieve Empirical Decision Times
//' 
//' This function only takes pLBA parameters and converts choice RTs to choice  
//' DTs in a matrix form. This is just a convenient function used interanlly by
//' \code{logLik_plba}. 
//'  
//' @param data a matrix with the first column stores choice (0==error; 
//' 1==correct), and the second column stores RTs in second.
//' @param pVec a vector storing pLBA model parameters. The sequence is 
//' critical. It is, A, b, muv1, muv2, t_ND, muw1, muw2, t_delay, sv, swt. 
//' See \code{\link{rplba}} for an explanation of the parameters.
//' @return A two-element list with decision time for accumulator 1 and 2.
//' @export
//' @examples
//' data(lba)
//' pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51,  
//'           muw2=3.69, t_delay=0.31,  sv=1, swt=0.5)
//'            
//' DTs <- cpda::choiceDT(data.matrix(d), pVec)
//'            
// [[Rcpp::export]]
arma::mat choiceDT(arma::mat data, arma::vec pVec) {
  // A  b muv1  muv2  t_ND  muw1  muw2 t_delay  sv  swt 
  // 0  1    2     3     4     5     6       7   8    9 
  // block_1 = c(0, 1, 2, 3, 4); ## A, b, muv1, muv2, t_ND
  // block_2 = c(5, 6, 7);       ## muw1, muw2, t_delay
  
  arma::vec choice = data.col(0) ;       // Response Time Block; 0==error; 
  arma::vec rt     = data.col(1) ;       // 1==correct
  arma::uvec idxX  = find(choice == 1) ; // error
  arma::uvec idxO  = find(choice == 2) ; // correct
  arma::vec timeX  = rt.rows(idxX) ;     // X RT; accumulator X
  arma::vec timeO  = rt.rows(idxO) ;     // O RT; accumulator O
  
  // This only place the sequence is matter.
  arma::vec DT1 = timeX - pVec[4] ; // X acc; substract off t0 
  arma::vec DT2 = timeO - pVec[4] ; // O acc
  
  arma::vec response1(DT1.n_elem);
  arma::vec response2(DT2.n_elem);
  response1.fill(1);
  response2.fill(2);
  
  arma::vec DT       = arma::join_cols(DT1, DT2);
  arma::vec response = arma::join_cols(response1, response2);
  arma::mat out      = arma::join_horiz(response, DT); 
  return out;
}


