#include <cpda.hpp>

// Protect against compilers without OpenMP; e.g., OS X  clang
#ifdef _OPENMP
#include <omp.h>
#endif

void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
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

//' @rdname rplba
//' @export
// [[Rcpp::export]]
arma::mat rplba1(int n, arma::vec pVec) {
  //  A  b  v1   v2  w1   w2  sv  rD  swt t0
  //  0  1   2    3   4    5   6   7    8  9 
  arma::vec x1(n), x2(n), v1(n), v2(n), dt1(n), dt2(n);  
  arma::vec choice(n), winDT(n), undone(n);
  double T0 = pVec[7] + pVec[8]; // rate delay + switch  

  // Stage 1 LBA
  for(int i=0; i<n; i++) 
  { 
    v1[i]  = rtn_scalar(pVec[2], pVec[6], 0, INFINITY) ; // acc1 (X)
    v2[i]  = rtn_scalar(pVec[3], pVec[6], 0, INFINITY) ; // acc2 (O)
    x1[i]  = Rf_runif(0, pVec[0]);
    x2[i]  = Rf_runif(0, pVec[0]);
    dt1[i] = (pVec[1] - x1[i]) / v1[i];
    dt2[i] = (pVec[1] - x2[i]) / v2[i];
    choice[i] = (dt1[i] < dt2[i]) ? 1 : 2;
    winDT[i]  = (dt1[i] < dt2[i]) ? dt1[i] : dt2[i];
    undone[i] = (winDT[i] <= T0) ? false : true;
  }
  
  int n2         = arma::accu(undone);
  arma::uvec idx = find(undone);
  arma::vec x1s2 = x1.elem(idx);  
  arma::vec x2s2 = x2.elem(idx);
  arma::vec v1s2 = v1.elem(idx);  
  arma::vec v2s2 = v2.elem(idx); 
  arma::vec dt1s2(n2), dt2s2(n2);

  //  A  b  v1   v2  w1   w2  sv  rD  swt t0
  //  0  1   2    3   4    5   6   7    8  9 
  // Stage 2 LBA
  for (int j=0; j<n2; j++) 
  {
    dt1s2[j] = (pVec[1] - (x1s2[j] + T0*v1s2[j])) / 
      (rtn_scalar(pVec[4], pVec[6], 0, INFINITY));
    dt2s2[j] = (pVec[1] - (x2s2[j] + T0*v2s2[j])) / 
      (rtn_scalar(pVec[5], pVec[6], 0, INFINITY));
    choice[idx[j]] = (dt1s2[j] < dt2s2[j]) ? 1 : 2;
    winDT[idx[j]]  = (dt1s2[j] < dt2s2[j]) ? (dt1s2[j] + T0) : (dt2s2[j] + T0);
  }
  
  arma::vec RT = winDT + pVec[9]; // Add t0 
  arma::mat out = arma::join_horiz(choice, RT);
  return out;
}

//' @rdname rplba
//' @export
// [[Rcpp::export]]
arma::mat rplba2(int n, arma::vec pVec) {
  /* A1  A2  b1  b2  v1  v2  w1  w2  sv1  sv2 sw1  sw2  rD swt  t0
      0   1   2   3   4   5   6   7    8    9  10   11  12  13  14 */
  arma::vec x1(n), x2(n), v1(n), v2(n), dt1(n), dt2(n);  
  arma::vec choice(n), winDT(n), undone(n);
  double T0 = pVec[12] + pVec[13]; // delay + switch
  
  // Stage 1 LBA
  for(int i=0; i<n; i++) 
  { 
    v1[i]  = rtn_scalar(pVec[4], pVec[8],  0, INFINITY) ; // acc1 (X)
    v2[i]  = rtn_scalar(pVec[5], pVec[9], 0, INFINITY) ; // acc2 (O)
    x1[i]  = Rf_runif(0, pVec[0]);
    x2[i]  = Rf_runif(0, pVec[1]);
    dt1[i] = (pVec[2] - x1[i]) / v1[i];
    dt2[i] = (pVec[3] - x2[i]) / v2[i];
    choice[i] = (dt1[i] < dt2[i]) ? 1 : 2;
    winDT[i]  = (dt1[i] < dt2[i]) ? dt1[i] : dt2[i];
    undone[i] = (winDT[i] <= T0) ? false : true;
  }
  
  int n2         = arma::accu(undone);
  arma::uvec idx = find(undone);
  arma::vec x1s2 = x1.elem(idx);  
  arma::vec x2s2 = x2.elem(idx);
  arma::vec v1s2 = v1.elem(idx);  
  arma::vec v2s2 = v2.elem(idx); 
  arma::vec dt1s2(n2), dt2s2(n2);
  
  /* A1  A2  b1  b2  v1  v2  w1  w2  sv1  sv2 sw1  sw2  rD swt  t0
      0   1   2   3   4   5   6   7    8    9  10   11  12  13  14 */
  
  // Stage 2 LBA
  for (int j=0; j<n2; j++) 
  {
    dt1s2[j] = (pVec[2] - (x1s2[j] + T0*v1s2[j])) / 
      (rtn_scalar(pVec[6], pVec[10], 0, INFINITY));
    dt2s2[j] = (pVec[3] - (x2s2[j] + T0*v2s2[j])) / 
      (rtn_scalar(pVec[7], pVec[11], 0, INFINITY));
    choice[idx[j]] = (dt1s2[j] < dt2s2[j]) ? 1 : 2;
    winDT[idx[j]]  = (dt1s2[j] < dt2s2[j]) ? (dt1s2[j] + T0) : (dt2s2[j] + T0);
  }
  
  arma::vec RT = winDT + pVec[14]; // Add t0 
  arma::mat out = arma::join_horiz(choice, RT);
  return out;
}

//' Generate Random Choice Response Times using pLBA Model
//'
//' This function uses two-accumulator piecewise LBA model to generate random 
//' choice RTs. There are 3 variants: \code{rplba}, \code{rplba1}, and 
//' \code{rplba2}. Each of them has a corresponding R version,
//' \code{rplbaR}, \code{rplbaR1}, and \code{rplbaR2}, for the purpose of 
//' speed testing. Because the difference of random number generators in C and 
//' R, they do not generate exactly identical RTs.  When generating large 
//' enough observations, the distributions generated by R and C will match.
//' 
//' The main function \code{rplba} implements a more flexible  
//' version of pLBA random number generator than the other two. It uses the  
//' following parameterisation (order matter):
//'  
//' \itemize{
//' \item \bold{\emph{A1}} accumulator 1 start-point upper bound. \code{A} is  
//' the upper bound of the interval \code{[0, A]}, which is used by an uniform 
//' distribution to generate a start-point. Average amount of 
//' prior evidence (i.e., before accumulation process even begins) across trials
//' is \code{A/2}.
//' \item \bold{\emph{A2}} accumulator 2 start-point upper bound.
//' \item \bold{\emph{B1}} accumulator 1 traveling distance. Note this is not 
//' a decision threshold!. LBA convention denotes decision threshold/caution as 
//' b (lowercase) and traveling distance as B (uppercase). \code{B=b-A} is 
//' the traveling distance, and \code{b-A/2} is a measure of average 
//' \emph{decision caution}.
//' \item \bold{\emph{B2}} accumulator 2 traveling distance.
//' \item \bold{\emph{C1}} the amount of traveling distance change for 
//' accumulator 1 at the stage 2. 
//' \item \bold{\emph{C2}} the amount of traveling distance change for 
//' accumulator 2 at the stage 2.
//' \item \bold{\emph{v1}} accumulator 1 drift rate, stage 1
//' \item \bold{\emph{v2}} accumulator 2 drift rate, stage 1
//' \item \bold{\emph{w1}} accumulator 1 drift rate, stage 2
//' \item \bold{\emph{w2}} accumulator 2 drift rate, stage 2
//' \item \bold{\emph{sv1}} accumulator 1 drift rate standard deviation, 
//' stage 1.
//' \item \bold{\emph{sv2}} accumulator 2 drift rate standard deviation, 
//' stage 1.
//' \item \bold{\emph{sw1}} accumulator 1 drift rate standard deviation, 
//' stage 2.
//' \item \bold{\emph{sw2}} accumulator 2 drift rate standard deviation, 
//' stage 2.
//' \item \bold{\emph{rD}} the delay duration while stage 1 drift rate switches
//' to stage 2 drift rate 
//' \item \bold{\emph{tD}} the delay duration while stage 1 threshold switches
//' to stage 2 threshold 
//' \item \bold{\emph{swt}} switch time, usually determined by experimental 
//' design. 
//' \item \bold{\emph{t0}} non-decision time in second. 
//' }
//' 
//' \code{rplba1} uses the following parameterisation: 
//' 
//' \itemize{
//' \item \bold{\emph{A}} a common start-point interval for both accumulators. 
//' \item \bold{\emph{b}} a common response threshold for both accumulators.  
//' \item \bold{\emph{v1}} accumulator 1 drift rate, stage 1
//' \item \bold{\emph{v2}} accumulator 2 drift rate, stage 1
//' \item \bold{\emph{w1}} accumulator 1 drift rate, stage 2
//' \item \bold{\emph{w2}} accumulator 2 drift rate, stage 2
//' \item \bold{\emph{sv}} a common standard deviation for both accumulators
//' \item \bold{\emph{rD}} a delay period while drift rate switch to a 
//' second stage process   
//' \item \bold{\emph{swt}} switch time, usually determined by experimental 
//' design
//' \item \bold{\emph{t0}} non-decision time in second. 
//' }
//' 
//' \code{rplba2} uses the following parameterisation: 
//' 
//' \itemize{
//' \item \bold{\emph{A1}} start-point interval of the accumulator 1. 
//' \item \bold{\emph{A2}} start-point interval of the accumulator 2. 
//' \item \bold{\emph{b1}} accumulator 1 response threshold. 
//' \item \bold{\emph{b2}} accumulator 2 response threshold. 
//' \item \bold{\emph{v1}} accumulator 1 drift rate, stage 1
//' \item \bold{\emph{v2}} accumulator 2 drift rate, stage 1
//' \item \bold{\emph{w1}} accumulator 1 drift rate, stage 2
//' \item \bold{\emph{w2}} accumulator 2 drift rate, stage 2
//' \item \bold{\emph{sv1}} the standard deviation of accumulator 1 drirt rate
//' during stage 1. 
//' \item \bold{\emph{sv2}} the standard deviation of accumulator 2 drirt rate
//' during stage 1. 
//' \item \bold{\emph{sw1}} the standard deviation of accumulator 1 drirt rate
//' during stage 2. 
//' \item \bold{\emph{sw2}} the standard deviation of accumulator 2 drirt rate
//' during stage 2. 
//' \item \bold{\emph{rD}} a delay period while drift rate switch to a 
//' second stage process   
//' \item \bold{\emph{swt}} switch time, usually determined by experimental 
//' design
//' \item \bold{\emph{t0}} non-decision time in second. 
//' }
//'
//' @param n number of observations. Must be an integer 
//' @param pVec a numeric vector storing pLBA model parameters. The sequence is
//' critical. See details for the sequence. 
//' @return A \code{n x 2} matrix with a first column storing choices and second 
//' column storing response times.
//' @export
//' @examples
//' ################
//' ## Example 1  ##
//' ################
//' pVec3.1 <- c(A1=1.51, A2=1.51, B1=1.2, B2=1.2, C1=.3, C2=.3, v1=3.32, 
//'              v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.1,
//'              tD=.1, swt=0.5, t0=0.08)
//' pVec3.2 <- c(A1=1.51, A2=1.51, B1=1.2, B2=1.2, C1=.3, C2=.3, v1=3.32, 
//'              v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.1,
//'              tD=.15, swt=0.5, t0=0.08)
//' pVec3.3 <- c(A1=1.51, A2=1.51, B1=1.2, B2=1.2, C1=.3, C2=.3, v1=3.32, 
//'              v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.15,
//'              tD=.1, swt=0.5, t0=0.08)
//' 
//' n <- 1e5
//' set.seed(123); system.time(dat5.1 <- cpda::rplbaR(n, pVec3.1))
//' set.seed(123); system.time(dat5.2 <- cpda::rplbaR(n, pVec3.2))
//' set.seed(123); system.time(dat5.3 <- cpda::rplbaR(n, pVec3.3))
//' set.seed(123); system.time(dat6.1 <- cpda::rplba( n, pVec3.1))
//' set.seed(123); system.time(dat6.2 <- cpda::rplba( n, pVec3.2))
//' set.seed(123); system.time(dat6.3 <- cpda::rplba( n, pVec3.3))
//' tmp5.1 <- data.frame(choice=factor(dat5.1[,1]), rt=dat5.1[,2])
//' tmp5.2 <- data.frame(choice=factor(dat5.2[,1]), rt=dat5.2[,2])
//' tmp5.3 <- data.frame(choice=factor(dat5.3[,1]), rt=dat5.3[,2])
//' tmp6.1 <- data.frame(choice=factor(dat6.1[,1]), rt=dat6.1[,2])
//' tmp6.2 <- data.frame(choice=factor(dat6.2[,1]), rt=dat6.2[,2])
//' tmp6.3 <- data.frame(choice=factor(dat6.3[,1]), rt=dat6.3[,2])
//'   
//' tmp5.1$fun <- "R"
//' tmp5.2$fun <- "R"
//' tmp5.3$fun <- "R"
//' tmp6.1$fun <- "C"
//' tmp6.2$fun <- "C"
//' tmp6.3$fun <- "C"
//' 
//' tmp5.1$vec <- "1"
//' tmp5.2$vec <- "2"
//' tmp5.3$vec <- "3"
//' tmp6.1$vec <- "1"
//' tmp6.2$vec <- "2"
//' tmp6.3$vec <- "3"
//' 
//' df <- rbind(tmp5.1, tmp5.2, tmp5.3, tmp6.1, tmp6.2, tmp6.3)
//' df$fun <- factor(df$fun)
//'         
//' ## Show R and C functions produce almost identical distributions        
//' \dontrun{ 
//' ## Set up a colour palette 
//' cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
//'         "#D55E00", "#CC79A7")
//'         
//' require(ggplot2)
//' ggplot(data=df, aes(x = rt, fill=fun, color=fun)) +
//'   geom_density(alpha=0.2) +
//'   facet_grid(vec~ choice) +
//'   scale_fill_manual(values=cb)
//'   
//' ## Or you can use lattice or base graphics  
//' require(lattice)
//' histogram( ~rt | vec+choice+fun, data=df, breaks="fd", type="density", 
//'            xlab="Response Time (s)", 
//'            panel=function(x, ...) {
//'                  panel.histogram(x, ...)
//'                  panel.densityplot(x, darg=list(kernel="gaussian"),...)
//'   })
//' }
//'              
//' par(mfrow=c(3,2))
//' hist(tmp5.1[tmp5.1$choice==1,"rt"], breaks="fd", col="gray", freq=FALSE, 
//'        xlab="RT (s)", main="pLBA-Choice 1")
//' lines(density(tmp6.1[tmp6.1$choice==1,"rt"]), col="red", lty="dashed",  lwd=1.5)
//'   
//' hist(tmp5.1[tmp5.1$choice==2,"rt"], breaks="fd", col="gray", freq=FALSE, 
//'        xlab="RT (s)", main="pLBA-Choice 2")
//' lines(density(tmp6.1[tmp6.1$choice==2,"rt"]), col="red", lty="dashed",  lwd=1.5)
//'   
//' #############
//' hist(tmp5.2[tmp5.2$choice==1,"rt"], breaks="fd", col="gray", freq=FALSE, 
//'        xlab="RT (s)", main="pLBA-Choice 1")
//' lines(density(tmp6.2[tmp6.2$choice==1,"rt"]), col="red", lty="dashed",  lwd=1.5)
//'     
//' hist(tmp5.2[tmp5.2$choice==2,"rt"], breaks="fd", col="gray", freq=FALSE, 
//'          xlab="RT (s)", main="pLBA-Choice 2")
//' lines(density(tmp6.2[tmp6.2$choice==2,"rt"]), col="red", lty="dashed",  lwd=1.5)
//'     
//' #############
//' hist(tmp5.3[tmp5.3$choice==1,"rt"], breaks="fd", col="gray", freq=FALSE, 
//'          xlab="RT (s)", main="pLBA-Choice 1")
//' lines(density(tmp6.3[tmp6.3$choice==1,"rt"]), col="red", lty="dashed",  lwd=1.5)
//'       
//' hist(tmp5.3[tmp5.3$choice==2,"rt"], breaks="fd", col="gray", freq=FALSE, 
//'            xlab="RT (s)", main="pLBA-Choice 2")
//' lines(density(tmp6.3[tmp6.3$choice==2,"rt"]), col="red", lty="dashed",  lwd=1.5)
//' par(mfrow=c(1,1))
//' 
//' ################
//' ## Example 2  ##
//' ################
//' pVec1 <- c(A=1.51, b=2.7, v1=3.32, v2=2.24,  w1=1.51,  w2=3.69,  
//'            sv=1, rD=0.31, swt=0.5, t0=0.08)
//'  
//' pVec2 <- c(A1=1.51, A2=1.51, b1=2.7, b2=2.7, v1=3.32, v2=2.24,
//'            w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, sw2=1, rD=0.31, 
//'            swt=0.5, t0=0.08)
//'            
//' system.time(dat1 <- cpda::rplba1( n, pVec1))
//' system.time(dat2 <- cpda::rplba2( n, pVec2))
//' system.time(dat3 <- cpda::rplbaR1(n, pVec1))
//' system.time(dat4 <- cpda::rplbaR2(n, pVec2))
//' 
//' tmp1 <- data.frame(choice=factor(dat1[,1]), rt=dat1[,2])
//' tmp2 <- data.frame(choice=factor(dat2[,1]), rt=dat2[,2])
//' tmp3 <- data.frame(choice=factor(dat3[,1]), rt=dat3[,2])
//' tmp4 <- data.frame(choice=factor(dat4[,1]), rt=dat4[,2])
//' tmp1$fun <- "rplba1"
//' tmp2$fun <- "rplba2"
//' tmp3$fun <- "rplba1-R"
//' tmp4$fun <- "rplba2-R"
//' tmp0 <- rbind(tmp1, tmp2, tmp3, tmp4)
//' tmp0$fun <- factor(tmp0$fun)
//'   
//' \dontrun{   
//' require(ggplot2)
//' ggplot(data = tmp0, aes(x = rt, fill=fun, color=fun)) +
//'     geom_density(alpha=0.2) +
//'     facet_grid(.~ choice) +
//'     scale_fill_manual(values=cb)
//' }
//'     
// [[Rcpp::export]]
arma::mat rplba(int n, arma::vec pVec) {
   // A1  A2  B1  B2  C1  C2  v1  v2  w1  w2  sv1  sv2 sw1  sw2  rD  tD  swt t0
   //  0   1   2   3   4   5   6   7   8   9   10   11  12   13  14  15  16  17 
   arma::vec x1(n), x2(n), v1(n), v2(n), dt1(n), dt2(n);  
   arma::vec choice(n), winDT(n), undone1(n);
   // Calculate two switch times
   double swt_r = pVec[14] + pVec[16]; // rate delay + switch
   double swt_b = pVec[15] + pVec[16]; // threshold delay + switch
   double b1 = pVec[0] + pVec[2];      // orignal threshold for acc 1
   double b2 = pVec[1] + pVec[3];      // orignal threshold for acc 2
   double c1 = b1 + pVec[4];           // changed threshold for acc 1
   double c2 = b2 + pVec[5];           // changed threshold for acc 2 
   double swt1, swt2;                  // swt mutation
   bool a0=false, a1=false, a2=false;
   
   // Determine which switching condition occurs, depending on rate delay and
   // threshold delay. Both parameters are feed from a sampler.
   if (swt_r == swt_b) {       // condition 0: rate and thresold change co-occur
     a0 = true; 
     swt1 = swt_r;
     swt2 = swt_r;
   } else if (swt_b < swt_r) { // condition 1: threshold change occurs early 
     a1 = true;
     swt1 = swt_b;
     swt2 = swt_r;
   } else { // (swt_b > swt_r); condition 2: rate chage occurs early 
     a2 = true; 
     swt1 = swt_r;
     swt2 = swt_b;
   }
   
   double swtD = swt2-swt1;

   for(int i=0; i<n; i++)           // Stage 1 LBA
   {
       v1[i]  = rtn_scalar(pVec[6], pVec[10], 0, INFINITY) ; // acc1 rate
       v2[i]  = rtn_scalar(pVec[7], pVec[11], 0, INFINITY) ; // acc2 rate
       x1[i]  = Rf_runif(0, pVec[0]);   
       x2[i]  = Rf_runif(0, pVec[1]);
       dt1[i] = (b1 - x1[i]) / v1[i]; // Race in original distances
       dt2[i] = (b2 - x2[i]) / v2[i];
       choice[i]  = (dt1[i] < dt2[i]) ? 1 : 2;
       winDT[i]   = (dt1[i] < dt2[i]) ? dt1[i] : dt2[i];
       undone1[i] = (winDT[i] <= swt1) ? false : true;
   }

  // A1  A2  B1  B2  C1  C2  v1  v2  w1  w2  sv1  sv2 sw1  sw2  rD  tD  swt t0
  //  0   1   2   3   4   5   6   7   8   9  10    11  12   13  14  15  16  17 
   int n2 = arma::accu(undone1);  
   arma::uvec idx = find(undone1);
   arma::vec x1s2 = x1.elem(idx);
   arma::vec x2s2 = x2.elem(idx);
   arma::vec v1s2 = v1.elem(idx); // acc1 observations with S1 drift rate not finished 
   arma::vec v2s2 = v2.elem(idx); // acc2 observations with S1 drift rate not finished 
   arma::vec dt1s2(n2), dt2s2(n2), B1s2(n2), B2s2(n2), w1(n2), w2(n2);
   arma::vec undone2(n2);

   for (int j=0; j<n2; j++)  // Stage 2 LBA; start time at swt1
   {
     // Calculate distance left to travel for those unfinished observations
     // threshold - (start-point + shorter swt * stage1 drift)  
     if (a2) { 
       // If rate change occurs earlier than threshold change 
       // (swt1==swt_r < swt_b), threshold stays at b.
       // swt1*v1s2[j] calculates distances already travelled, so v1s2, 
       // instead of w1, B2s2 ditto.
       B1s2[j] = b1 - (x1s2[j] + swt1*v1s2[j]); 
       B2s2[j] = b2 - (x2s2[j] + swt1*v2s2[j]); 
     } else { 
       // If swt_b shorter or equal to swt_r, threshold has already changed.
       // swt1==swt_b < swt_r.
       // Total travel distance changes to c. That is, b1 to c1 and b2 to c2; 
       B1s2[j] = c1 - (x1s2[j] + swt1*v1s2[j]); 
       B2s2[j] = c2 - (x2s2[j] + swt1*v2s2[j]); 
     }
     
     if (a1) {  
       // If threshold change occurs earlier than rate change (swt_b < swt_r)
       // continue to use stage 1 drift rates, because rates stay the same.
       w1[j] = v1s2[j];
       w2[j] = v2s2[j];
     } else {
       // swt_b is longer than or equal to swt_r. rate change has occurred.
       w1[j] = rtn_scalar(pVec[8], pVec[12], 0, INFINITY);
       w2[j] = rtn_scalar(pVec[9], pVec[13], 0, INFINITY);
     }
     
     dt1s2[j]       = B1s2[j] / w1[j]; 
     dt2s2[j]       = B2s2[j] / w2[j]; 
     choice[idx[j]] = (dt1s2[j] < dt2s2[j]) ? 1 : 2;
     winDT[idx[j]]  = (dt1s2[j] < dt2s2[j]) ? (dt1s2[j]+swt1) : (dt2s2[j]+swt1);
     undone2[j]     = (winDT[idx[j]] < swt2) ? false : true;
   }
   
   // Stage 3 occurs only when swt_r and swt_b are unequal, namely a0=false
   if (!a0) {
     int n3 = arma::accu(undone2);  
     if(n3 > 0) {
       arma::uvec idx2 = find(undone2);
       arma::vec B1s3  = B1s2.elem(idx2);
       arma::vec B2s3  = B2s2.elem(idx2);
       arma::vec w1s3  = w1.elem(idx2); // acc1 observations with S2 drift rate not finished 
       arma::vec w2s3  = w2.elem(idx2); // acc2 observations with S2 drift rate not finished 
       arma::vec dt1s3(n3), dt2s3(n3), vvv1(n3), vvv2(n3), B1last(n3), B2last(n3);
       
       // Distance left to travel for those not finished
       // Distance left at end of stage 1 - further travel
       for (int k=0; k<n3; k++) { // Stage 3 LBA
         B1last[k] = B1s3[k] - swtD * w1s3[k];
         B2last[k] = B2s3[k] - swtD * w2s3[k];

         if (a1) {  // (swt_b < swt_r)
           // If threshold change occurs earlier than rate change, now
           // is the time for rates to change.
           vvv1[k] = rtn_scalar(pVec[8], pVec[12], 0, INFINITY);
           vvv2[k] = rtn_scalar(pVec[9], pVec[13], 0, INFINITY);
         } else {
           // If rate change occurs earlier than threshold change, 
           // inherit drift rates from previous stage and change threshold by
           // add C onto travel distances
           vvv1[k] = w1s3[k];
           vvv2[k] = w2s3[k];
           B1last[k] = B1last[k] + pVec[4];
           B2last[k] = B2last[k] + pVec[5];
         }
         
         dt1s3[k] = B1last[k] / vvv1[k];
         dt2s3[k] = B2last[k] / vvv2[k];
         choice[idx[idx2[k]]] = (dt1s3[k] < dt2s3[k]) ? 1 : 2;
         winDT[idx[idx2[k]]]  = (dt1s3[k] < dt2s3[k]) ? (dt1s3[k]+swt2) : (dt2s3[k]+swt2);
       }  
     }
   }
   
   arma::vec RT = winDT + pVec[17]; // Add t0 
   arma::mat out = arma::join_horiz(choice, RT);
   return out;
}
