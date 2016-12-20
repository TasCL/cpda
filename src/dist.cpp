#include <cpda.hpp>

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
//' @param ns number of simulations. 
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
double logLik_norm(arma::vec object, arma::vec pVec, int ns, 
                   double h=0, double m=0.8, double p=10) {
  // pVec[0] is mean, pVec[1] is sigma
  arma::vec y = arma::sort(object) ;
  // if (h==0) { h = bwNRD0(y, m); 
  // } else { h = m*h; } 
  // int nd         = y.n_elem ;
  arma::vec yhat = pVec[0]+pVec[1] * arma::randn(ns) ; // simulation
  double LL      = logLik_fft(y, yhat, h, m, p) ;
  return LL ;
}

//' @rdname logLik_norm
//' @export
// [[Rcpp::export]]
Rcpp::List logLik_norm2(arma::vec object, arma::vec pVec, int ns,
                        double h=0, double m=0.8, double p=10) {
  arma::vec y = arma::sort(object) ;
  arma::vec yhat = pVec[0]+pVec[1] * arma::randn(ns) ; // simulation
  Rcpp::List LL  = logLik_fft2(y, yhat, h, m, p) ;
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
//' critical. It is, A, b, muv1, muv2, t_ND, muw1, muw2, t_delay, sv, swt. 
//' @param ns number of simulations. Usually \code{cpda} can handle up to 1e6. 
//' Use \code{gpda} if large simulation is required.
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
//' pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51, 
//'           muw2=3.69, t_delay=0.31, sv=1, swt=0.5)
//'           
//' data(lba)
//' dMat <- data.matrix(d); head(dMat)
//' 
//' set.seed(123)
//' samp   <- cpda::rplba(1e5, pVec); head(samp)
//' h   <- 0.8*bw.nrd0(samp[,2]); h
//' set.seed(123)
//' ll0 <- cpda::logLik_plba(dMat, pVec, 1e5, h=h, m=1); ll0
//' ## [1] -7254.58
//' set.seed(123)
//' llList <- cpda::logLik_plba2(dMat, pVec, 1e5, h=h, m=1); 
//' str(llList)
//' 
//' ## List of 4
//' ## $ LL      : num -7255
//' ## $ PDF     : num [1:1000, 1] -23 -23 -23 -23 -23 ...
//' ## $ z       : num [1:2048, 1] 0.175 0.179 0.184 0.188 0.192 ...
//' ## $ PDF_hist: num [1:2048, 1] 0 0 0 0 0 0 0 0 0 0 ...
//' 
//' ## Here I demonstrate use rplba and logLik_fft to produce identical 
//' ## result
//' set.seed(123)
//' samp   <- cpda::rplba(1e5, pVec); head(samp)
//' DTMat  <- cpda::choiceDT(dMat, pVec)
//' time1  <- sort(DTMat[DTMat[,1] == 1, 2])
//' time2  <- sort(DTMat[DTMat[,1] == 2, 2])
//' time1_ <- sort(samp[samp[,1] == 1, 2])
//' time2_ <- sort(samp[samp[,1] == 2, 2])
//' ll1 <- cpda::logLik_fft(time1, time1_, h=h, m=1); ll1 ## [1] -5492.192
//' ll2 <- cpda::logLik_fft(time2, time2_, h=h, m=1); ll2 ## [1] -1762.388
//' print(ll1+ll2)
//' ## [1] -7254.58
//'          
// [[Rcpp::export]]
double logLik_plba(arma::mat object, arma::vec pVec, int ns,
                   double h=0, double m=0.8, double p=10) {
  // A  b muv1  muv2  t_ND  muw1  muw2 t_delay  sv  swt
  // 0  1    2     3     4     5     6       7   8    9  
  // block_1 = c(0, 1, 2, 3, 4); ## A, b, muv1, muv2, t_ND
  // block_2 = c(5, 6, 7);       ## muw1, muw2, t_delay

  arma::mat data   = choiceDT(object, pVec) ;
  arma::vec R      = data.col(0);
  arma::vec dt     = data.col(1);
  arma::vec time1  = arma::sort(dt.rows(find(R == 1))); // acc1 data DT; 
  arma::vec time2  = arma::sort(dt.rows(find(R == 2))); // acc2 data DT; 
  
  arma::mat sim    = rplba(ns, pVec); 
  arma::vec R_     = sim.col(0); // choice; 1 or 2
  arma::vec dt_    = sim.col(1); // rplba returns DTs
  arma::vec time1_ = arma::sort(dt_.rows(find(R_ == 1))); // acc1 sim DT; 
  arma::vec time2_ = arma::sort(dt_.rows(find(R_ == 2))); // acc2 sim DT; 

  double LL1 = logLik_fft(time1, time1_, h, m, p) ;
  double LL2 = logLik_fft(time2, time2_, h, m, p) ;

  return (LL1+LL2);
}

//' @rdname logLik_plba
//' @export
// [[Rcpp::export]]
Rcpp::List logLik_plba2(arma::mat object, arma::vec pVec, int ns,
                   double h=0, double m=0.8, double p=10) {
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
  
  Rcpp::List out1 = logLik_fft2(time1, time1_, h, m, p) ;
  Rcpp::List out2 = logLik_fft2(time2, time2_, h, m, p) ;
  double LL1       = out1["LL"];
  double LL2       = out2["LL"];
  arma::vec PDF1   = out1["PDF"];
  arma::vec PDF2   = out2["PDF"];
  arma::vec z1     = out1["z"];
  arma::vec z2     = out2["z"];
  arma::vec PDF_hist1 = out1["PDF_hist"];
  arma::vec PDF_hist2 = out2["PDF_hist"];
  
  double LL          = LL1 + LL2;
  arma::vec PDF      = arma::join_cols(PDF1, PDF2);
  arma::vec z        = arma::join_cols(z1, z2);
  arma::vec PDF_hist = arma::join_cols(PDF_hist1, PDF_hist2);
  
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
//' @param ns number of simulations.  
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
double logLik_lba(arma::mat object, arma::vec pVec, int ns, double h=0, 
                  double m=0.8, double p=10) {
  arma::vec R   = object.col(0);
  arma::vec rt  = object.col(1);
  arma::vec dt1 = arma::sort(rt.rows(find(R == 1)) - pVec[8]); // acc1 data RT; 
  arma::vec dt2 = arma::sort(rt.rows(find(R == 2)) - pVec[9]); // acc1 data RT; 

  arma::mat sim    = rlba(ns, pVec); // rlba returns rts
  arma::vec R_     = sim.col(0);     // response; 1 vs 2
  arma::vec rt_    = sim.col(1);     // DTs
  arma::vec dt1_   = arma::sort(rt_.rows(find(R_ == 1)) - pVec[8]); // substract t01
  arma::vec dt2_   = arma::sort(rt_.rows(find(R_ == 2)) - pVec[9]); // substract t02

  double LL1 = logLik_fft(dt1, dt1_, h, m, p) ; // return summed, logged lik
  double LL2 = logLik_fft(dt2, dt2_, h, m, p) ; // return summed, logged lik
  
  return (LL1+LL2);
}

//' @rdname logLik_lba
//' @export
// [[Rcpp::export]]
Rcpp::List logLik_lba2(arma::mat object, arma::vec pVec, int ns,
                  double h=0, double m=0.8, double p=10) {
  arma::vec R   = object.col(0);
  arma::vec rt  = object.col(1);
  arma::vec dt1 = arma::sort(rt.rows(find(R == 1)) - pVec[8]); // acc1 data RT; 
  arma::vec dt2 = arma::sort(rt.rows(find(R == 2)) - pVec[9]); // acc1 data RT; 
  
  arma::mat sim    = rlba(ns, pVec); // rlba returns rts
  arma::vec R_     = sim.col(0);     // response; 1 vs 2
  arma::vec rt_    = sim.col(1);     // DTs
  arma::vec dt1_   = arma::sort(rt_.rows(find(R_ == 1)) - pVec[8]); // substract t01
  arma::vec dt2_   = arma::sort(rt_.rows(find(R_ == 2)) - pVec[9]); // substract t02

  Rcpp::List out1  = logLik_fft2(dt1, dt1_, h, m, p) ; 
  Rcpp::List out2  = logLik_fft2(dt2, dt2_, h, m, p) ; 
  double LL1       = out1["LL"];
  double LL2       = out2["LL"];
  arma::vec PDF1   = out1["PDF"];
  arma::vec PDF2   = out2["PDF"];
  arma::vec z1     = out1["z"];
  arma::vec z2     = out2["z"];
  arma::vec PDF_hist1 = out1["PDF_hist"];
  arma::vec PDF_hist2 = out2["PDF_hist"];
  
  double LL          = LL1 + LL2;
  arma::vec PDF      = arma::join_cols(PDF1, PDF2);
  arma::vec z        = arma::join_cols(z1, z2);
  arma::vec PDF_hist = arma::join_cols(PDF_hist1, PDF_hist2);
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("LL")        = LL,
    Rcpp::Named("PDF")       = PDF,
    Rcpp::Named("z")         = z,
    Rcpp::Named("PDF_hist")  = PDF_hist);

  return (out);
}
