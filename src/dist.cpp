#include <cpda.hpp>

//' Compute FFT Log-likelihood for a Gaussian Distribution
//'
//' This is a wrapper function to approximate a Gaussian distribution, using
//' KDE-FFT method. To retrieve more outputs and replicate Holmes's example 1,
//' please use \code{logLik_norm2}.
//'
//' @param object a vector storing empirical data
//' @param pVec parameter vector storing mean and standard deviation
//' @param ns number of simulations. 
//' @return Log-likelihood; plus PDF, grid centers, and PDF_hist
//' @export
//' @examples
//' pVec <- c(mu=5, sigma=1)
//' y    <- sort(rnorm(1e5, pVec[1], pVec[2]))
//'
//' ll1 <- logLik_norm(y, pVec, 1e5)
//' ll2 <- logLik_norm2(y, pVec, 1e5)
//' str(ll1); str(ll2)
//'
//' plot(ll2$z, ll2$PDF_hist, type="l", lty=2,
//' main="Normal Distribution",xlab="x",ylab="L(x|beta)")
//' lines(y, ll2$PDF, col="red", lwd = 1.5)
//'
// [[Rcpp::export]]
double logLik_norm(arma::vec object, arma::vec pVec, int ns) {
  // pVec[0] is mean, pVec[1] is sigma
  arma::vec y = arma::sort(object) ;
  int nd      = y.n_elem ;
  arma::vec yhat = pVec[0]+pVec[1] * arma::randn(ns) ; // simulation
  double LL      = logLik_fft(y, yhat) ;
  return LL ;
}

//' @rdname logLik_norm
//' @export
// [[Rcpp::export]]
Rcpp::List logLik_norm2(arma::vec object, arma::vec pVec, int ns) {
  arma::vec y = arma::sort(object) ;
  int nd      = y.n_elem ;
  arma::vec yhat = pVec[0]+pVec[1] * arma::randn(ns) ; // simulation
  
  Rcpp::List LL  = logLik_fft2(y, yhat) ;
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("LL")        = LL["LL"],
    Rcpp::Named("PDF")       = LL["PDF"],
    Rcpp::Named("z")         = LL["z"], // grid centers
    Rcpp::Named("PDF_hist")  = LL["PDF_hist"]) ;
  return out ;
}

//' Compute FFT Log-likelihood for a pLBA Model
//'
//' This is a wrapper function to approximate a pLBA distribution, using
//' KDE-FFT method. To retrieve more outputs, please use \code{logLik_pLBA2}.
//'
//' @param object a matrix storing empirical choice RT data. First column must
//' stores choices and second column stores RTs in seconds.
//' @param pVec a vector storing pLBA model parameters. The sequence is 
//' critical. It is, A, b, muv1, muv2, t_ND, muw1, muw2, t_delay, sv, swt. 
//' @param ns number of simulations. Usually \code{cpda} can handle up to 1e6. 
//' Use \code{gpda} if large simulation is required.
//' @return a summed, logged likelihood across trials and accumulators.
//' @export
//' @examples  
//' rm(list=ls())
//' data(lba)
//' dMat <- data.matrix(d)
//' head(dMat)
//' 
//' pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51, muw2=3.69,
//'           t_delay=0.31, sv=1, swt=0.5)
//' ## setting <- c(bandwidth=.02, ns=1e5, nmc=30, nchain=24, rp=.001, burnin=10,
//' ## nthin=3, start=1, gammaMult=2.38, report=100)
//' tmp0 <- cpda::logLik_plba(dMat, pVec, 1e5); tmp0
//'  
// [[Rcpp::export]]
double logLik_plba(arma::mat object, arma::vec pVec, int ns) {
  // A  b muv1  muv2  t_ND  muw1  muw2 t_delay  sv  swt
  // 0  1    2     3     4     5     6       7   8    9  
  // block_1 = c(0, 1, 2, 3, 4); ## A, b, muv1, muv2, t_ND
  // block_2 = c(5, 6, 7);       ## muw1, muw2, t_delay

  arma::mat data   = choiceDT(object, pVec) ;
  arma::vec R      = data.col(0);
  arma::vec dt     = data.col(1);
  arma::vec time1  = arma::sort(dt.rows(find(R == 1))); // acc1 data DT; 
  arma::vec time2  = arma::sort(dt.rows(find(R == 2))); // acc2 data DT; 
  
  arma::mat sim    = rplba(ns, pVec); // the only place pVec matter
  arma::vec R_     = sim.col(0); // response; 1 vs 2
  arma::vec dt_    = sim.col(1); // DTs // rplba returns dts
  arma::vec time1_ = arma::sort(dt_.rows(find(R_ == 1))); // acc1 sim DT; 
  arma::vec time2_ = arma::sort(dt_.rows(find(R_ == 2))); // acc2 sim DT; 

  double LL1 = logLik_fft(time1, time1_) ;
  double LL2 = logLik_fft(time2, time2_) ;

  return (LL1+LL2);
}

//' Compute FFT Log-likelihood for a LBA Model
//'
//' This is a wrapper function to approximate a LBA distribution, using
//' KDE-FFT method. To retrieve more outputs, please use \code{logLik_LBA2}.
//'
//' @param object a matrix storing empirical choice RT data. First column must
//' stores choices and second column stores RTs in seconds.
//' @param pVec pLBA parameter vector. The sequence is 
//' critical. It is, b1, b2, A1, A2, mu1, mu2, sigma1, sigma2, t01, and t02. 
//' @param ns number of simulations. Usually \code{cpda} can handle up to 1e6. 
//' Use \code{gpda} if large simulation is required.
//' @return a summed, logged likelihood across trials and accumulators.
//' @export
//' @examples  
//' rm(list=ls())
//' data(lba)
//' dMat <- data.matrix(d)
//' head(dMat)
//'  
//' pVec <- c(b1=1, b2=1, A1=.5, A2=.5, mu1=2.4, mu2=1.6, sigma1=1, sigma2=1.2,
//' t01=.5, t02=.5)
//' tmp0 <- cpda::logLik_lba(dMat, pVec, 1e5); tmp0
//'  
// [[Rcpp::export]]
double logLik_lba(arma::mat object, arma::vec pVec, int ns) {
  arma::mat data   = choiceDT(object, pVec) ;
  arma::vec R      = data.col(0);
  arma::vec dt     = data.col(1);
  arma::vec time1  = arma::sort(dt.rows(find(R == 1))); // acc1 data DT; 
  arma::vec time2  = arma::sort(dt.rows(find(R == 2))); // acc2 data DT; 
  
  arma::mat sim    = rlba(ns, pVec);  // rlba returns rts
  arma::vec R_     = sim.col(0); // response; 1 vs 2
  arma::vec rt_    = sim.col(1); // DTs
  arma::vec dt1_   = rt_.rows(find(R_ == 1)) - pVec[8]; // substract t01
  arma::vec dt2_   = rt_.rows(find(R_ == 2)) - pVec[9]; // substract t02
  arma::vec time1_ = arma::sort(dt1_); // acc1 sim DT; 
  arma::vec time2_ = arma::sort(dt2_); // acc2 sim DT; 
  
  double LL1 = logLik_fft(time1, time1_) ; // return summed, logged lik
  double LL2 = logLik_fft(time2, time2_) ; // return summed, logged lik
  
  return (LL1+LL2);
}

