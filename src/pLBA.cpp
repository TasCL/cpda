#include <cpda.hpp>

extern "C" void logLik_pw(double *y_, double *yhat_, int *ny, int *ns, 
                          double *h_, double *out);

arma::vec getEdges(arma::vec z) {
  // Get bin_edges
  double dt   = z[1] - z[0] ;
  arma::vec term1 = z - dt/2 ;
  arma::vec term2(1) ;
  double last = z[z.size()-1] + dt/2;
  term2.fill(last) ;
  arma::vec out = arma::join_cols(term1, term2) ;
  return out ;
}

arma::vec getFilter(double m, double M, double h) {
  // cannoical Gaussian kernel
  double tmp0    = 2 * arma::datum::pi * (std::pow(2, 10) / (M-m)) * 0.5 ;
  arma::vec tmp1 = arma::linspace<arma::vec>(0, 1, 1 + (std::pow(2, 10)/2)) ;
  arma::vec freq = tmp0 * tmp1 ;
  arma::vec s2   = arma::pow(freq, 2) ; // s^2 on p17
  double h2      = std::pow(h, 2) ;
  arma::vec fil0 = arma::exp(-0.5 * h2 * s2) ;
  arma::vec fil1 = arma::flipud(fil0.rows(1, (fil0.size() - 2)));
  arma::vec out  = arma::join_cols(fil0, fil1) ;
  return out ;
}

arma::vec density(arma::vec y, arma::vec be, double dt) {
  // y is yhat; be is binEdges;
  arma::uvec hc       = arma::histc(y, be) ;
  arma::vec bincount  = arma::conv_to<arma::vec>::from(hc);
  int ns              = arma::sum(bincount);
  arma::vec PDF_hist  = bincount / (dt * ns);
  arma::vec out       = PDF_hist.rows(0, (PDF_hist.size() - 2)) ;
  return out ;
}

//' Generate Random Choice Response Times Using LBA Model
//'
//' This function uses two-accumulator LBA model to generate random choice RTs.
//' 
//' \itemize{
//' \item \bold{\emph{b1}} response threshold for accumulator 1. 
//' \item \bold{\emph{b2}} response threshold for accumulator 2. 
//' \item \bold{\emph{A1}} starting point interval for accumulator 1. 
//' \item \bold{\emph{A2}} starting point interval for accumulator 2. 
//' \item \bold{\emph{mu1}} drift rate for accumulator 1.
//' \item \bold{\emph{mu2}} drift rate for accumulator 2. 
//' \item \bold{\emph{sigma1}} drift rate standard deviation for accumulator 1.
//' \item \bold{\emph{sigma2}} drift rate standard deviation for accumulator 2.
//' \item \bold{\emph{t01}} non-decision time for accumulator 1. 
//' \item \bold{\emph{t02}} non-decision time for accumulator 2. 
//' }
//'
//' @param n number of observations. Must be an integer
//' @param pVec a two-accumulator LBA parameter vector. The sequence is 
//' critical. It is, b1, b2, A1, A2, mu1, mu2, sigma1, sigma2, t01, and t02.
//' @return A matrix. First column stores choices and second column stores RTs 
//' in seconds
//' @export
//' @examples
//' n    <- 10
//' pVec <- c(b1=1, b2=1, A1=.5, A2=.5, mu1=2.4, mu2=1.6, sigma1=1, sigma2=1.2,
//'           t01=.5, t02=.5)
//' dat1 <- cpda::rlba(n, pVec)
//' str(dat1)
//' summary(dat1)
//' 
//' ## Compare to rtdists
//' ## dat2 <- rtdists::rLBA(n, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), 
//' ## sd_v=c(1,1.2))
//' ## str(dat2)
//' ## summary(dat2)
//' 
// [[Rcpp::export]]
arma::mat rlba(int n, arma::vec pVec) {
  arma::vec rts(n);
  arma::vec responses(n);
  //  pVec: b1 b2 A1 A2 mu1 mu2 sigma1 sigma2 t01 t02 
  //         0  1  2  3   4   5      6      7   8   9
  arma::vec b1     = arma::repmat(pVec.row(0), n, 1);
  arma::vec b2     = arma::repmat(pVec.row(1), n, 1);
  arma::vec A1     = arma::repmat(pVec.row(2), n, 1);
  arma::vec A2     = arma::repmat(pVec.row(3), n, 1);
  arma::vec mu1    = arma::repmat(pVec.row(4), n, 1);
  arma::vec mu2    = arma::repmat(pVec.row(5), n, 1);
  arma::vec sigma1 = arma::repmat(pVec.row(6), n, 1);
  arma::vec sigma2 = arma::repmat(pVec.row(7), n, 1);

  arma::vec x01 = A1 % arma::randu(n);
  arma::vec x02 = A2 % arma::randu(n);
  arma::vec v1  = mu1 + sigma1 % arma::randn(n);
  arma::vec v2  = mu2 + sigma2 % arma::randn(n);

  for (int i = 0; i < n; i++)
  {
    if(v1[i] < 0) {v1[i] = 1e6;}
    if(v2[i] < 0) {v2[i] = 1e6;}
  }

  arma::vec dt1 = (b1 - x01) / v1;
  arma::vec dt2 = (b2 - x02) / v2;

  for (int i = 0; i < n; i++)
  {
    if(dt1[i] < dt2[i]) {
      rts[i] = dt1[i] + pVec[8];
      responses[i] = 1;
    } else {
      rts[i] = dt2[i] + pVec[9];
      responses[i] = 2;
    }

  }
  
  arma::mat out = arma::join_horiz(responses, rts); 
  return out;
}

//' Generate Random Choice Decision Times
//'
//' Uses two-accumulator pLBA model to generate random choice DTs (not RTs!).
//' 
//' \itemize{
//' \item \bold{\emph{A}} starting point interval. A starting point is with the 
//' interval \code{[0, A]}. Average amount of prior evidence (ie before accumulation 
//' process even begins) across trials is \code{A/2}. 
//' \item \bold{\emph{b}} response threshold. \code{b-A/2} is a measure of 
//' \emph{response caution}.
//' \item \bold{\emph{muv1}} accumulator 1 drift rate, piece 1
//' \item \bold{\emph{muv2}} accumulator 2 drift rate, piece 1
//' \item \bold{\emph{t_ND}} non-decision time in second. 
//' \item \bold{\emph{muw1}} accumulator 1 drift rate, piece 2
//' \item \bold{\emph{muw2}} accumulator 2 drift rate, piece 2
//' \item \bold{\emph{t_delay}} delay time for the two-piece process   
//' \item \bold{\emph{sv}} a common standard deviation for all drift rates 
//' (muv1, muw1, muv2, muw2).
//' \item \bold{\emph{swt}} switch time, usually determined by experimental 
//' design).
//' }
//'
//' @param n number of observations. Must be an integer 
//' @param pVec a double vector storing pLBA model parameters. The sequence is 
//' critical. It is, A, b, muv1, muv2, t_ND, muw1, muw2, t_delay, sv, swt. 
//' @return A two-element list with estimated decision time for accumulator 1,
//' accumultor 2, empirical decision time for accumualtor 1 and accumulator 2.
//' @export
//' @examples
//' pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51, 
//'           muw2=3.69, t_delay=0.31, sv=1, swt=0.5)
//' ## setting <- c(bandwidth=.02, ns=1e5, nmc=30, nchain=24, rp=.001, burnin=10,
//' ## nthin=3, start=1, gammaMult=2.38, report=100)
//' DTs <- cpda::rplba(n=1e4, pVec)
//' summary(DTs)
//' head(DTs)
//' 
// [[Rcpp::export]]
arma::mat rplba(int n, arma::vec pVec) {
  /* A  b muv1  muv2  t_ND  muw1  muw2 t_delay  sv  swt; c(0, 1, 2, 3, 4)
     0  1    2     3     4     5     6       7   8    9; c(5, 6, 7)
     h ns nmc nchain rp burnin nthin start gammaMult report
     0  1   2      3  4      5     6     7         8      9 */
  double T0    = pVec[9] + pVec[7];  // switch time + the delay time 
  arma::vec x1 = pVec[0]*arma::randu(n); // X acc; Sample starting points
  arma::vec x2 = pVec[0]*arma::randu(n); // O acc 
  
  arma::vec v1(n); arma::vec w1(n);   // Sample drift rates
  arma::vec v2(n); arma::vec w2(n);
  // If n > 1e6, use omp. Maybe CPU can まわしげり GPU 
  for(int i=0; i<n; i++) // No negative drift rates
  { 
    v1[i] = rtn_scalar(pVec[1], pVec[8], 0, INFINITY) ; // X acc
    w1[i] = rtn_scalar(pVec[5], pVec[8], 0, INFINITY) ; // X acc
    v2[i] = rtn_scalar(pVec[3], pVec[8], 0, INFINITY) ; // O acc
    w2[i] = rtn_scalar(pVec[6], pVec[8], 0, INFINITY) ; // O acc
  }

  // Compute the time at which the two accumulators (v1==X & v2==O)
  // will terminate. v1 & v2 are two accumulators of the 1st piece
  // w1 & w2 are two accumulators of 2nd piece
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
  
  arma::mat out = arma::join_horiz(response, DT); 
  return out;
}

//' Retrieve Empirical Decision Times
//' 
//' This function takes pLBA parameters and converts choice RTs, stored in a 
//' matrix to choice DTs, stored in a list with two numeric vectors. 
//'  
//' @param data a double matrix with the first column stores choice (0==error; 
//' 1==correct), and the second column stores RTs in second.
//' @param pVec a double vector storing pLBA model parameters. The sequence is 
//' critical. It is, A, b, muv1, muv2, t_ND, muw1, muw2, t_delay, sv, swt. 
//' @return A two-element list with decision time for accumulator 1 and 2.
//' @export
//' @examples
//' data(lba)
//' pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51,  
//'           muw2=3.69, t_delay=0.31,  sv=1, swt=0.5)
//'            
//' DTs <- cpda::choiceDT(data.matrix(d), pVec)
//' 
//' ## Wrong sequence
//' pVecWrong <- c(A=1.51, muv1=3.32, muw1=1.51, muv2=2.24, muw2=3.69,
//'               t_delay=0.31, t_ND=0.08, b=2.7, sv=1, swt=0.5)
//' DTX <- cpda::choiceDT(data.matrix(d), pVecWrong)
//' head(DTX)
//'            
// [[Rcpp::export]]
arma::mat choiceDT(arma::mat data, arma::vec pVec) {
  // A  b muv1  muv2  t_ND  muw1  muw2 t_delay  sv  swt 
  // 0  1    2     3     4     5     6       7   8    9 
  // block_1 = c(0, 1, 2, 3, 4); ## A, b, muv1, muv2, t_ND
  // block_2 = c(5, 6, 7);       ## muw1, muw2, t_delay
  
  arma::vec choice = data.col(0) ;       // Response Time Block; 0==error; 
  arma::vec rt     = data.col(1) ;       // 1==correct
  arma::uvec idxX  = find(choice == 0) ; // error
  arma::uvec idxO  = find(choice == 1) ; // correct
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

//' Compute Log-likelihood Using KDE-based Fast Fourier Transform
//'
//' This function implements Holmes's (2015) KDE-FFT method to calculate
//' approximated probability density. The precision is set at 10 (ie 2^10).
//' \code{logLik_fft2} returns three additional vectors, PDF (discrete
//' probability densities for interpolated empirical data), grid centers, and
//' PDF_hist (discrete probability densities for the simulated data).
//'
//' @param y a vector storing empirical data (e.g., RTs)
//' @param yhat a vector storing simulated data (e.g., simualted RTs, using a
//' LBA model).
//' @return Log-likelihood
//' @references Holmes, W. (2015). A practical guide to the Probability Density
//' Approximation (PDA) with improved implementation and error characterization.
//' \emph{Journal of Mathematical Psychology}, \bold{68-69}, 13--24,
//' doi: http://dx.doi.org/10.1016/j.jmp.2015.08.006.
//' @export
//' @examples
//' ## Use piecewise LBA data as an example
//' data(lba)
//' logLik_fft(plba$DT1, plba$eDT1)
//' logLik_fft(plba$DT2, plba$eDT2)
//' tmp1 <- logLik_fft2(plba$DT1, plba$eDT1)
//' tmp2 <- logLik_fft2(plba$DT2, plba$eDT2)
//'
//' str(tmp1)
//' ## List of 4
//' ## $ LL      : num 278
//' ## $ PDF     : num [1:695, 1] 2.441 0.632 1.966 2.359 2.214 ...
//' ## $ z       : num [1:1024, 1] 0.0889 0.0903 0.0916 0.093 0.0943 ...
//' ## $ PDF_hist: num [1:1024, 1] 0 0 0 0 0 0 0 0 0 0 ...
//' str(tmp2)
//' ## List of 4
//' ## $ LL      : num 42.9
//' ## $ PDF     : num [1:305, 1] 1.879 0.965 1.834 0.326 0.921 ...
//' ## $ z       : num [1:1024, 1] 0.067 0.0696 0.0722 0.0748 0.0774 ...
//' ## $ PDF_hist: num [1:1024, 1] 0 0 0 0 0 0 0 0 0 0 ...
// [[Rcpp::export]]
double logLik_fft(arma::vec y, arma::vec yhat) {
  int nd = y.n_elem ;
  int ns = yhat.n_elem ;
  double h = bwNRD0(y, 0.8);

  double z0 = std::min(y.min(), yhat.min()) - 3 * h;
  double z1 = std::max(y.max(), yhat.max()) + 3 * h;
  arma::vec z        = arma::linspace<arma::vec>(z0, z1, std::pow(2, 10)) ;
  arma::vec binEdges = getEdges(z) ;
  arma::vec filter   = getFilter(z0, z1, h) ;  // Gauss filter
  double dt          = z[1] - z[0] ;
  arma::vec PDF_tmp ;

  arma::vec signal   = density(yhat, binEdges, dt) ;
  arma::cx_vec PDF_fft      = arma::fft(signal) ;
  arma::cx_vec PDF_fft_filt = filter % PDF_fft ;
  arma::vec PDF_smoothed = arma::real(arma::ifft(PDF_fft_filt)) ;
  arma::interp1(z, PDF_smoothed, y, PDF_tmp);
  arma::vec PDF = pmax(PDF_tmp, std::pow(10, -5)) ;
  /* int nSignal        = signal.n_elem;
     int nFilter        = filter.n_elem;
     if (nSignal != nFilter) 
          fprintf(stderr, "signal & filter have unequal sizes\n");  */
  double LL = arma::accu(arma::log(PDF));
  return LL ;
}

//' @rdname logLik_fft
//' @export
// [[Rcpp::export]]
Rcpp::List logLik_fft2(arma::vec y, arma::vec yhat) {
  int nd   = y.n_elem ;
  int ns   = yhat.n_elem ;
  //double h = nd==1 ? h_ : bwNRD0(y, 0.8);
  double h = bwNRD0(y, 0.8);

  double z0   = std::min(y.min(), yhat.min()) - 3 * h;
  double z1   = std::max(y.max(), yhat.max()) + 3 * h;
  arma::vec z = arma::linspace<arma::vec>(z0, z1, std::pow(2, 10)) ;
  arma::vec binEdges = getEdges(z) ;
  arma::vec filter   = getFilter(z0, z1, h) ;  // Gauss filter
  double dt          = z[1] - z[0] ;
  arma::vec signal   = density(yhat, binEdges, dt) ;
  /* int nSignal        = signal.n_elem; // 2^10
  int nFilter        = filter.n_elem; // 2^10
  if (nSignal != nFilter) fprintf(stderr, "signal & filter have unequal sizes\n");
  */
  arma::cx_vec PDF_fft      = arma::fft(signal) ;
  arma::cx_vec PDF_fft_filt = filter % PDF_fft ;
  arma::vec PDF_smoothed    = arma::real(arma::ifft(PDF_fft_filt)) ;

  arma::vec PDF_tmp;   // Interpolate the grid likelihood to the data
  arma::interp1(z, PDF_smoothed, y, PDF_tmp);
  arma::vec PDF = pmax(PDF_tmp, std::pow(10, -5)) ;
  double LL     = arma::accu(arma::log(PDF));

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("LL")        = LL,
    Rcpp::Named("PDF")       = PDF,
    Rcpp::Named("z")         = z,
    Rcpp::Named("PDF_hist")  = signal);
  return out ;
}

void logLik_pw(double *y_, double *yhat_, int *ny, int *ns, double *h_, 
               double *out) {
  arma::vec yhat = getVec(yhat_ ,ns);
  for(int i=0; i<*ny; i++)
  {
    out[i] = gaussian(y_[i], yhat, *h_) ;
  }
}

