arma::mat rplba_wh(int n, arma::vec pVec) {
  /* A  b muv1  muv2  t_ND  muw1  muw2 t_delay  sv  swt; c(0, 1, 2, 3, 4)
  0  1    2     3     4     5     6       7   8    9; c(5, 6, 7) */
  // double T0    = pVec[9] + pVec[7];      // switch time + the delay time 
  double T0    = pVec[9];  // + pVec[7];    // switch time + the delay time 
  arma::vec x1 = pVec[0]*arma::randu(n); // X acc; Sample starting points
  arma::vec x2 = pVec[0]*arma::randu(n); // O acc 
  arma::vec v1(n); arma::vec w1(n);      // Drift rate containers
  arma::vec v2(n); arma::vec w2(n);
  
  for(int i=0; i<n; i++) // Use my rtnorm functions to exclude negative v's 
  { 
    v1[i] = rtn_scalar(pVec[2], pVec[8], 0, INFINITY) ; // X acc
    w1[i] = rtn_scalar(pVec[5], pVec[8], 0, INFINITY) ; // X acc
    v2[i] = rtn_scalar(pVec[3], pVec[8], 0, INFINITY) ; // O acc
    w2[i] = rtn_scalar(pVec[6], pVec[8], 0, INFINITY) ; // O acc
  }
  
  // Compute the time at which the two accumulators (v1==X & v2==O)
  // will terminate. v1 & v2 are the rates of two accumulators during 1st stage 
  // w1 & w2 are the rates of two accumulators during 2nd stage
  arma::vec A1P1 = (pVec[1]-x1) / v1; // X acc; accumulator 1, piece 1
  arma::vec A2P1 = (pVec[1]-x2) / v2; // O acc; accumulator 2, piece 1
  arma::vec A1P2 = T0 + (pVec[1]-x1-v1*T0) / w1; // X acc; accumulator 1, piece 2
  arma::vec A2P2 = T0 + (pVec[1]-x2-v2*T0) / w2; // O acc; accumulator 2, piece 2
  
  // pre switch (< T0) termination
  arma::uvec A1WinB = find((A1P1 < T0) && (A1P1 < A2P1)) ; // A1 wins before
  arma::uvec A2WinB = find((A2P1 < T0) && (A2P1 < A1P1)) ; // A2 wins before
  
  // post switch (> T0) termination
  arma::uvec A1WinA = find((A1P1 > T0) && (A2P1 > T0) && (A1P2 < A2P2) &&
    (A1P2 < 1e4)) ; // A1 wins after switch
  arma::uvec A2WinA = find((A1P1 > T0) && (A2P1 > T0) && (A2P2 < A1P2) &&
    (A2P2 < 1e4)) ; // A2 wins after switch
  
  arma::vec DT1 = arma::join_cols(A1P1.elem( A1WinB ), A1P2.elem( A1WinA ));
  arma::vec DT2 = arma::join_cols(A2P1.elem( A2WinB ), A2P2.elem( A2WinA ));
  
  arma::vec response1(DT1.n_elem);
  arma::vec response2(DT2.n_elem);
  response1.fill(1);
  response2.fill(2);
  
  arma::vec DT       = arma::join_cols(DT1, DT2);
  arma::vec response = arma::join_cols(response1, response2);
  // arma::mat out = arma::join_horiz(response, DT); 
  
  arma::vec DT_ = DT + pVec[7];
  arma::mat out = arma::join_horiz(response, DT_);
  
  return out;
}
