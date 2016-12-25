#include <cpda.hpp>

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



//' Compute Summed Log-likelihood for a Gaussian Distribution
//'
//' This is a wrapper function to approximate a Gaussian distribution, using
//' KDE-FFT method. To retrieve more outputs, use \code{logLik_norm2}.
//'
//' \code{logLik_norm2} is identical to \code{logLik_norm}, except returning 
//' four elements:
//' 
//' \itemize{
//' \item \bold{\emph{LL}}, summed, logged likelihood. This is the same as the 
//' return value from \code{logLik_plba}. 
//' \item \bold{\emph{PDF}}, a numeric vector storing \emph{logged} 
//' probability densities for individual data point.
//' \item \bold{\emph{z}}, a numeric vector storing centre points of the 
//' simulated histogram (i.e., grid centre)
//' \item \bold{\emph{PDF_hist}} a numeric vector stoing the count of simulated
//' data point in each bin 
//' } 
//' @param object a vector storing empirical data
//' @param pVec parameter vector storing mean and standard deviation
//' @param ns number of simulations. Default is 1e5.
//' @param h KDE bandwidth. If not input been enter, the default is 
//' Sliverman's rule of thumb; otherwise the function uses the input from
//' the user; 
//' @param m a multiplier to adjust proportationally h. Default is 0.8
//' @param p a precision parameter defines the number of grid as power of 2.
//' Default value is 10 (i.e., 2^10).
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
//' main="Normal Distribution",xlab="x",ylab="L(x|theta)")
//' lines(y, exp(ll2$PDF), col="red", lwd = 1.5)
//'
// [[Rcpp::export]]
double logLik_norm(arma::vec object, arma::vec pVec, int ns=1e5, 
                   double h=0, double m=0.8, double p=10) {
  // pVec[0] is mean, pVec[1] is sigma
  arma::vec yhat = pVec[0]+pVec[1] * arma::randn(ns) ; // simulation
  double LL      = logLik_fft(object, yhat, h, m, p) ;
  return LL ;
}

//' @rdname logLik_norm
//' @export
// [[Rcpp::export]]
Rcpp::List logLik_norm2(arma::vec object, arma::vec pVec, int ns=1e5,
                        double h=0, double m=0.8, double p=10) {
  arma::vec yhat = pVec[0]+pVec[1] * arma::randn(ns) ; // simulation
  Rcpp::List LL  = logLik_fft2(object, yhat, h, m, p) ;
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("LL")        = LL["LL"],
    Rcpp::Named("PDF")       = LL["PDF"],
    Rcpp::Named("z")         = LL["z"], // grid centers
    Rcpp::Named("PDF_hist")  = LL["PDF_hist"]) ;
  return out ;
}

//' Compute Summed Log-likelihood for a pLBA Model
//'
//' This is a wrapper function to approximate a pLBA distribution, using
//' KDE-FFT method. To retrieve more outputs, use \code{logLik_plba2}.
//'
//' \code{logLik_plba2} is identical to \code{logLik_plba}, except returning 
//' four elements:
//' 
//' \itemize{
//' \item \bold{\emph{LL}}, summed, logged likelihood. This is the same as the 
//' return value from \code{logLik_plba}. 
//' \item \bold{\emph{PDF}}, a numeric vector storing \emph{logged}  
//' probability densities for individual data point.
//' \item \bold{\emph{z}}, a numeric vector storing centre points of the 
//' simulated histogram (i.e., grid centre)
//' \item \bold{\emph{PDF_hist}} a numeric vector stoing the count of simulated
//' data point in each bin 
//' } 
//' 
//' @param object a matrix storing empirical choice RT data. First column must
//' stores choices and second column stores RTs in seconds.
//' @param pVec a vector storing pLBA model parameters. The sequence is 
//' critical. it is A1, A2, b1, b2, v1, v2, w1, w2, sv1, sv2, sw1, sw2,
//' rD, swt, and t0. 
//' @param ns number of simulations. Default is 1e5.
//' @param h KDE bandwidth. If not input been enter, the default is 
//' Sliverman's rule of thumb; otherwise the function uses the input from
//' the user; 
//' @param m a multiplier to adjust proportationally h. Default is 0.8
//' @param p a precision parameter defines the number of grid as power of 2.
//' Default value is 10 (i.e., 2^10).
//' @return a summed, logged likelihood across trials and accumulators. 
//' \code{logLik_plba2} returns three more elements.
//' @export
//' @examples  
//' ###################
//' ## Example 1     ##
//' ###################
//' ## I demonstrate how to use the build-in logLik_plba to calculate
//' ## pLBA densities. set.seed(123) is to produce the same result.
//' rm(list=ls())
//' pVec <- c(A1=1.51, A2=1.51, b1=2.7, b2=2.7, v1=3.32, v2=2.24,
//'             w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.31, 
//'             swt=0.5, t0=0.08)
//'   
//' data(lba)
//' d$Response <- ifelse(d$Response==0, 1, 2)
//' dMat <- data.matrix(d); head(dMat)
//' n <- 1e5
//'   
//' set.seed(123)
//' ll0 <- cpda::logLik_plba(dMat, pVec, 1e5, h=h, m=1); ll0
//'       
//' set.seed(123)
//' llList <- cpda::logLik_plba2(dMat, pVec, 1e5, h=h, m=1); str(llList)
//' ## List of 4
//' ## $ LL      : num 327
//' ## $ PDF     : num [1:1000, 1:4] 1 1 1 1 1 1 1 1 1 1 ...
//' ## $ z       : num [1:2048, 1] 0.13 0.135 0.14 0.146 0.151 ...
//' ## $ PDF_hist: num [1:2048, 1] 0 0 0 0 0 ...
//'         
//' ###################
//' ## Example 2     ##
//' ###################
//' ## Secondly, I demonstrate use rplba2 and logLik_fft to produce identical 
//' ## result
//' rm(list=ls())
//' x <- cbind(rep(1:2,each=100),rep(seq(.5,2,length.out=100),2))
//' pVec <- c(A1=1.51, A2=1.51, b1=2.7, b2=2.7, v1=3.32, v2=2.24,
//'             w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.31, 
//'             swt=0.5, t0=0.08)
//' dt1 <- x[x[,1]==1,2] - pVec[15]
//' dt2 <- x[x[,1]==2,2] - pVec[15]
//' 
//' set.seed(123)
//' n <- 1e5
//' samp <- cpda::rplba2(n, pVec); head(samp)
//' dt1_ <- samp[samp[,1]==1,2] - pVec[15]
//' dt2_ <- samp[samp[,1]==2,2] - pVec[15]
//' logLik_fft(dt1, dt1_) + logLik_fft(dt2, dt2_)
//'   
//' set.seed(123)
//' cpda::logLik_plba(x, pVec)
//'   
//' fft1 <- cpda::logLik_plba2(x, pVec)
//' str(fft1)
//' 
//' ## choice RT PDF
//' tmp0 <- fft1$PDF[fft1$PDF[,1] == 1, 2]
//' tmp1 <- sort(x[x[,1]==1,2])
//' all(tmp0==tmp1)
//'   
//' tmp3 <- fft1$PDF[fft1$PDF[,1] == 2,2]
//' tmp4 <- sort(x[x[,1]==1,2])
//' all(tmp3==tmp4)
//'               
// [[Rcpp::export]]
double logLik_plba(arma::mat object, arma::vec pVec, int ns=1e5,
                   double h=0, double m=0.8, double p=10) {
  /* A1  A2  b1  b2  v1  v2  w1  w2  sv1  sv2 sw1  sw2  rD swt  t0
      0   1   2   3   4   5   6   7    8    9  10   11  12  13  14 */
  arma::vec choice = object.col(0) ;  // 1==error; 2==correct 
  arma::vec rt     = object.col(1) ;        
  arma::vec dt1    = rt.rows(find(choice == 1)) - pVec[14]; // acc1 
  arma::vec dt2    = rt.rows(find(choice == 2)) - pVec[14]; // acc2
  arma::mat sim    = rplba2(ns, pVec); // choice RT
  arma::vec choice_= sim.col(0); // choice; 1 or 2
  arma::vec rt_    = sim.col(1); // rplba2 returns RTs
  arma::vec dt1_   = rt_.rows(find(choice_ == 1)) - pVec[14]; // acc1 sim DT; 
  arma::vec dt2_   = rt_.rows(find(choice_ == 2)) - pVec[14]; // acc2 sim DT; 
  return (logLik_fft(dt1, dt1_, h, m, p) + logLik_fft(dt2, dt2_, h, m, p));
}

//' @rdname logLik_plba
//' @export
// [[Rcpp::export]]
Rcpp::List logLik_plba2(arma::mat object, arma::vec pVec, int ns=1e5,
                   double h=0, double m=0.8, double p=10) {
  /* A1  A2  b1  b2  v1  v2  w1  w2  sv1  sv2 sw1  sw2  rD swt  t0
   0   1   2   3   4   5   6   7    8    9  10   11  12  13  14 */
  arma::vec choice = object.col(0) ;  // 1==error; 2==correct 
  arma::vec rt     = object.col(1) ;        
  arma::vec dt1    = rt.rows(find(choice == 1)) - pVec[14]; // acc1 
  arma::vec dt2    = rt.rows(find(choice == 2)) - pVec[14]; // acc2
  arma::mat sim    = rplba2(ns, pVec); // choice RT
  arma::vec choice_= sim.col(0); // choice; 1 or 2
  arma::vec rt_    = sim.col(1); // rplba2 returns RTs
  arma::vec dt1_   = rt_.rows(find(choice_ == 1)) - pVec[14]; // acc1 sim DT; 
  arma::vec dt2_   = rt_.rows(find(choice_ == 2)) - pVec[14]; // acc2 sim DT; 
  Rcpp::List ll1   = logLik_fft2(dt1, dt1_, h, m, p) ;
  Rcpp::List ll2   = logLik_fft2(dt2, dt2_, h, m, p) ;
  
  double LL1       = ll1["LL"];
  double LL2       = ll2["LL"];
  arma::mat PDF1   = ll1["PDF"]; // dt, PDF
  arma::mat PDF2   = ll2["PDF"]; 
  arma::vec z1     = ll1["z"];
  arma::vec z2     = ll2["z"];
  arma::vec hist1  = ll1["PDF_hist"];
  arma::vec hist2  = ll2["PDF_hist"];
  
  arma::vec choice1(PDF1.n_rows); choice1.fill(1); 
  arma::vec choice2(PDF2.n_rows); choice2.fill(2);
  arma::vec rtOut1  = PDF1.col(0) + pVec[14];
  arma::vec rtOut2  = PDF2.col(0) + pVec[14];
  arma::mat choiceRT1 = arma::join_horiz(choice1, rtOut1); 
  arma::mat choiceRT2 = arma::join_horiz(choice2, rtOut2);  
  arma::mat PDFMat1 = arma::join_horiz(choiceRT1, PDF1.col(1)); // choice RT, PDF
  arma::mat PDFMat2 = arma::join_horiz(choiceRT2, PDF2.col(1));

  double LL          = LL1 + LL2;
  arma::mat PDF      = arma::join_cols(PDFMat1, PDFMat2);
  arma::mat z        = arma::join_cols(z1, z2);
  arma::mat PDF_hist = arma::join_cols(hist1, hist2);

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("LL")        = LL,
    Rcpp::Named("PDF")       = PDF,
    Rcpp::Named("z")         = z,
    Rcpp::Named("PDF_hist")  = PDF_hist);
  return (out);
}

//' Compute Summed Log-likelihood for a LBA Model
//'
//' This is a wrapper function to approximate the densities of a LBA model, 
//' using \code{logLik_fft}. To retrieve more outputs, use \code{logLik_lba2}.
//' 
//' \code{logLik_lba2} is identical to \code{logLik_lba}, except the former 
//' returns three more elements: 
//'  
//' \itemize{
//' \item \bold{\emph{LL}}, summed, logged likelihood. This is the same as the 
//' return value from \code{logLik_lba}. 
//' \item \bold{\emph{PDF}}, a numeric vector storing logged probability 
//' densities for individual data point.
//' \item \bold{\emph{z}}, a numeric vector storing centre points of the 
//' simulated histogram (i.e., grid centre)
//' \item \bold{\emph{PDF_hist}} a numeric vector stoing the count of simulated
//' data point in each bin 
//' } 
//'
//' @param object a matrix storing empirical choice RT data. First column must
//' be choices (1 & 2) and second column must bve RTs in seconds.
//' @param pVec LBA parameter vector. The sequence is critical. It is b1, b2, 
//' A1, A2, mu1, mu2, sigma1, sigma2, t01, and t02. 1 and 2 stand for 
//' accumulator 1 and 2, respectively.  
//' @param ns number of simulations. Default is 1e5.
//' @param h KDE bandwidth. If not given, the default \code{h} is 
//' Sliverman's rule of thumb; otherwise the function uses the input from
//' the user; 
//' @param m a multiplier to adjust proportationally h. Default is 0.8
//' @param p a precision parameter defines the number of grid as power of 2.
//' Default value is 10 (i.e., 2^10).
//' @return a summed, logged likelihood across trials and accumulators. 
//' \code{logLik_plba2} returns three more elements.
//' @export
//' @examples  
//' ## Demo identical outputs from manually building rlba and logLik_fft, 
//' ## comparing to logLik_lba 
//' 
//' ## First retrieve example data set
//' data(lba)  
//' d$R <- ifelse(d$Response==0, 1, 2)  ## convert 0 & 1 accumulator to 1 & 2
//' dMat <- data.matrix(data.frame(R=d$R, RT=d$ResponseTime))
//' head(dMat)
//'  
//' ## LBA parameter vector. The sequence is critical. 
//' pVec <- c(b1=1, b2=1, A1=.5, A2=.5, mu1=2.4, mu2=1.6, sigma1=1, sigma2=1.2,
//'           t01=.5, t02=.5)
//'            
//' set.seed(123)  ## make sure using identical simulations
//' samp <- cpda::rlba(1e5, pVec); head(samp)
//' h    <- 0.8*bw.nrd0(samp[,2]); h
//' 
//' ## logLik_lba simulates internally, so set.seed ahead to make sure using
//' ## identical simulations
//' set.seed(123)  
//' tmp0 <- cpda::logLik_lba(dMat, pVec, 1e5, h, 1); tmp0 ## -3496.88
//' 
//' ## Manually calculate empirical DTs for accumulator 1 and accumualtor 2
//' dt1  <- sort(dMat[dMat[,1] == 1, 2]) - pVec[9]
//' dt2  <- sort(dMat[dMat[,1] == 2, 2]) - pVec[10]
//' dt1_ <- sort(samp[samp[,1] == 1, 2]) - pVec[9]
//' dt2_ <- sort(samp[samp[,1] == 2, 2]) - pVec[10]
//'     
//' ## Calculating log-likelihoods separately for each accumulator
//' ll1 <- cpda::logLik_fft(dt1, dt1_, h, 1); ll1
//' ll2 <- cpda::logLik_fft(dt2, dt2_, h, 1); ll2
//' print(ll1+ll2) ## -3496.88
//'  
// [[Rcpp::export]]
double logLik_lba(arma::mat object, arma::vec pVec, int ns=1e5, double h=0, 
                  double m=0.8, double p=10) {
  arma::vec R   = object.col(0);
  arma::vec rt  = object.col(1);
  arma::vec dt1 = rt.rows(find(R == 1)) - pVec[8]; // acc1 data RT; 
  arma::vec dt2 = rt.rows(find(R == 2)) - pVec[9]; // acc1 data RT; 

  arma::mat sim    = rlba(ns, pVec); // rlba returns rts
  arma::vec R_     = sim.col(0);     // response; 1 vs 2
  arma::vec rt_    = sim.col(1);     // DTs
  arma::vec dt1_   = rt_.rows(find(R_ == 1)) - pVec[8]; // substract t01
  arma::vec dt2_   = rt_.rows(find(R_ == 2)) - pVec[9]; // substract t02

  double LL1 = logLik_fft(dt1, dt1_, h, m, p) ; // return summed, logged lik
  double LL2 = logLik_fft(dt2, dt2_, h, m, p) ; // return summed, logged lik
  
  return (LL1+LL2);
}

//' @rdname logLik_lba
//' @export
// [[Rcpp::export]]
Rcpp::List logLik_lba2(arma::mat object, arma::vec pVec, int ns=1e5,
                  double h=0, double m=0.8, double p=10) {
  arma::vec R   = object.col(0);
  arma::vec rt  = object.col(1);
  arma::vec dt1 = rt.rows(find(R == 1)) - pVec[8]; // acc1 data RT; 
  arma::vec dt2 = rt.rows(find(R == 2)) - pVec[9]; // acc1 data RT; 
  arma::mat sim    = rlba(ns, pVec); // rlba returns rts
  arma::vec R_     = sim.col(0);     // response; 1 vs 2
  arma::vec rt_    = sim.col(1);     // DTs
  arma::vec dt1_   = rt_.rows(find(R_ == 1)) - pVec[8]; // substract t01
  arma::vec dt2_   = rt_.rows(find(R_ == 2)) - pVec[9]; // substract t02

  Rcpp::List ll1  = logLik_fft2(dt1, dt1_, h, m, p) ; 
  Rcpp::List ll2  = logLik_fft2(dt2, dt2_, h, m, p) ; 
  
  double LL1       = ll1["LL"];
  double LL2       = ll2["LL"];
  arma::mat PDF1   = ll1["PDF"]; // dt, PDF
  arma::mat PDF2   = ll2["PDF"]; 
  arma::vec z1     = ll1["z"];
  arma::vec z2     = ll2["z"];
  arma::vec hist1  = ll1["PDF_hist"];
  arma::vec hist2  = ll2["PDF_hist"];
  
  arma::vec choice1(PDF1.n_rows); choice1.fill(1); 
  arma::vec choice2(PDF2.n_rows); choice2.fill(2);
  arma::vec rtOut1  = PDF1.col(0) + pVec[8];
  arma::vec rtOut2  = PDF2.col(0) + pVec[9];
  arma::mat choiceRT1 = arma::join_horiz(choice1, rtOut1); 
  arma::mat choiceRT2 = arma::join_horiz(choice2, rtOut2);  
  arma::mat PDFMat1 = arma::join_horiz(choiceRT1, PDF1.col(1)); // choice RT, PDF
  arma::mat PDFMat2 = arma::join_horiz(choiceRT2, PDF2.col(1));
  
  double LL          = LL1 + LL2;
  arma::mat PDF      = arma::join_cols(PDFMat1, PDFMat2);
  arma::mat z        = arma::join_cols(z1, z2);
  arma::mat PDF_hist = arma::join_cols(hist1, hist2);
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("LL")        = LL,
    Rcpp::Named("PDF")       = PDF,
    Rcpp::Named("z")         = z,
    Rcpp::Named("PDF_hist")  = PDF_hist);
  return (out);
}
