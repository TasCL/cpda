#include <cpda.hpp>

//' @export
// [[Rcpp::export]]
double logLik_lba(arma::mat y, arma::vec pVec, int n) {
  arma::vec RT  = y.col(0);
  arma::vec RT1 = RT.rows(find(y.col(1) == 1)); // y.col(1) == R
  arma::vec RT2 = RT.rows(find(y.col(1) == 2)); 
  double out;

  // When pVec is abnormal, there would be chances that rlba generates sRT1, but
  // 0 sRT2. This is the occassion that when the function blows up itself.
  arma::mat sim  = rlba(n, pVec);
  arma::vec sRT  = sim.col(0);
  arma::vec sRT1 = sRT.rows(find(sim.col(1) == 1)); // sim.col(1) == R
  arma::vec sRT2 = sRT.rows(find(sim.col(1) == 2)); 
  
  arma::uvec sRT1_idx = arma::find(sRT1 < 4);
  arma::uvec sRT2_idx = arma::find(sRT2 < 4);
  arma::vec sRT1_trimmed = sRT1.elem(sRT1_idx);
  arma::vec sRT2_trimmed = sRT2.elem(sRT2_idx);

  if (sRT1.n_elem == 0) { 
    out = logLik_fft(RT2, sRT2_trimmed, 0, 0.8, 10, n);
  } else if (sRT2.n_elem == 0) {
    out = logLik_fft(RT1, sRT1_trimmed, 0, 0.8, 10, n);
  } else {
    out = logLik_fft(RT1, sRT1_trimmed, 0, 0.8, 10, n) + 
      logLik_fft(RT2, sRT2_trimmed, 0, 0.8, 10, n);
  }
  return (out);
}

