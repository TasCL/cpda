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

//' Generate Random Choice Decision Times using pLBA Model
//'
//' This function use two-accumulator piecewise LBA model to generate random 
//' choice DTs (not RTs!). A parallel version \code{rplba_omp} enables parallel 
//' random-number generation (used when n > 1e6).
//' 
//' \itemize{
//' \item \bold{\emph{A}} start point interval. \code{A} is the upper bound 
//' of the interval \code{[0, A]}. An uniform distribution with bound, 
//' \code{[0, A]} is used to generated a start point. Average amount of prior 
//' evidence (ie before accumulation process even begins) across trials is 
//' \code{A/2}. 
//' \item \bold{\emph{b}} response threshold. \code{b-A/2} is a measure of 
//' average \emph{response caution}.
//' \item \bold{\emph{muv1}} accumulator 1 drift rate, stage 1
//' \item \bold{\emph{muv2}} accumulator 2 drift rate, stage 1
//' \item \bold{\emph{t_ND}} non-decision time in second. 
//' \item \bold{\emph{muw1}} accumulator 1 drift rate, stage 2
//' \item \bold{\emph{muw2}} accumulator 2 drift rate, stage 2
//' \item \bold{\emph{t_delay}} delay time for the two-stage process   
//' \item \bold{\emph{sv}} a common standard deviation for all drift rates 
//' (muv1, muw1, muv2, muw2).
//' \item \bold{\emph{swt}} switch time, usually determined by experimental 
//' design.
//' }
//'
//' @param n number of observations. Must be an integer 
//' @param pVec a double vector storing pLBA model parameters. The sequence is 
//' critical. It is, A, b, muv1, muv2, t_ND, muw1, muw2, t_delay, sv, swt. 
//' @return A n x 2 matrix with a first column storing choices and second 
//' column storing DTs.
//' @export
//' @examples
//' ## pVec stands for parameter vector
//' pVec <- c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, muw1=1.51, 
//'           muw2=3.69, t_delay=0.31, sv=1, swt=0.5)
//' n <- 10
//' 
//' ## A R-rplba function
//' \dontrun{ 
//' rplba <- function(n, p=c(A=1.51, b=2.7, muv1=3.32, muv2=2.24, t_ND=0.08, 
//'                         muw1=1.51, muw2=3.69, t_delay=0.31, sv=1, swt=0.5)) 
//' {
//' require(msm)
//' ## Stage 1 LBA
//' v1 <- msm::rtnorm(n,p["muv1"],p["sv"],0); 
//' v2 <- msm::rtnorm(n,p["muv2"],p["sv"],0); 
//' sp <- matrix(runif(2*n,0,p["A"]),nrow=2)
//'     
//' ## Race
//' dt <- rbind((p["b"]-sp[1,])/v1,(p["b"]-sp[2,])/v2)
//' 
//' ## dt[dt<0] <- Inf
//' choice <- apply(dt,2,which.min)
//' rt <- dt[cbind(choice,1:n)]
//'     
//' ## Which are finished?
//' done <- rt <= p["swt"]
//' n2 <- sum(!done)
//'       
//' ## Distance left to travel for those not finished
//' B1 <- p["b"] - (sp[1,!done] + p["swt"]*v1[!done]) 
//' B2 <- p["b"] - (sp[2,!done] + p["swt"]*v2[!done]) 
//'         
//' ## Stage 2 LBA
//' w1 <- msm::rtnorm(n2,p["muw1"],p["sv"],0)
//' w2 <- msm::rtnorm(n2,p["muw2"],p["sv"],0)
//'           
//' ## Race  
//' dt <- rbind(B1/w1,B2/w2)
//' ## dt[dt<0] <- Inf
//' choice[!done] <- apply(dt,2,which.min)
//' rt[!done] <- p["swt"]+dt[cbind(choice[!done],1:n2)]
//'           
//' ## save results
//' cbind(choice=choice,rt=p["t_delay"]+rt)
//'             
//' }
//' }
//' 
//' system.time(DT1 <- rplba(n=1e2, p=pVec))
//' head(DT1)
//' 
//' \dontrun{
//' system.time(DT2 <- cpda::rplba(n=1e2, pVec=pVec))
//' head(DT2)
//' }
//' 
//' \dontrun{
//' ## Convert to data frame for plotting
//' tmp1 <- data.frame(choice=factor(DT1[,1]), dt=DT1[,2])
//' tmp2 <- data.frame(choice=factor(DT2[,1]), dt=DT2[,2])
//' summary(tmp1[tmp1$choice==1,])
//' summary(tmp2[tmp2$choice==1,])
//' 
//' summary(tmp1[tmp1$choice==2,])
//' summary(tmp2[tmp2$choice==2,])
//' 
//' tmp1$package <- "cpda"
//' tmp2$package <- "R"
//' tmp0 <- rbind(tmp1, tmp2)
//' tmp0$package <- factor(tmp0$package)
//' 
//' require(lattice)
//' lattice::histogram( ~ dt | package + choice, data = tmp0, aspect = 1,
//' xlab = "Decision Time (s)")
//'  
//' lattice::densityplot( ~ dt | package + choice, data = tmp0,   
//'          xlab = "Decisoin Time (s)", bw = stats::bw.nrd0(tmp0$dt))
//'}
// [[Rcpp::export]]
arma::mat rplba(int n, arma::vec pVec) {
  /* A  b muv1  muv2  t_ND  muw1  muw2 t_delay  sv  swt; c(0, 1, 2, 3, 4)
  0  1    2     3     4     5     6       7   8    9; c(5, 6, 7) */
  arma::vec x1(n), x2(n), v1(n), v2(n), dt1(n), dt2(n);  
  arma::vec choice(n), winDT(n), undone(n);
  
  // Stage 1 LBA
  for(int i=0; i<n; i++) 
  { 
    v1[i]  = rtn_scalar(pVec[2], pVec[8], 0, INFINITY) ; // acc1 (X)
    v2[i]  = rtn_scalar(pVec[3], pVec[8], 0, INFINITY) ; // acc2 (O)
    x1[i]  = Rf_runif(0, pVec[0]);
    x2[i]  = Rf_runif(0, pVec[0]);
    dt1[i] = (pVec[1] - x1[i]) / v1[i];
    dt2[i] = (pVec[1] - x2[i]) / v2[i];
    choice[i] = (dt1[i] < dt2[i]) ? 1 : 2;
    winDT[i]  = (dt1[i] < dt2[i]) ? dt1[i] : dt2[i];
    undone[i] = (winDT[i] <= pVec[9]) ? false : true;
  }
  
  int n2         = arma::accu(undone);
  arma::uvec idx = find(undone);
  arma::vec x1s2 = x1.elem(idx);  arma::vec x2s2 = x2.elem(idx);
  arma::vec v1s2 = v1.elem(idx);  arma::vec v2s2 = v2.elem(idx); 
  arma::vec dt1s2(n2), dt2s2(n2);
  
  // Stage 2 LBA
  for (int j=0; j<n2; j++) 
  {
    dt1s2[j] = (pVec[1] - (x1s2[j] + pVec[9]*v1s2[j])) / 
      (rtn_scalar(pVec[5], pVec[8], 0, INFINITY));
    dt2s2[j] = (pVec[1] - (x2s2[j] + pVec[9]*v2s2[j])) / 
      (rtn_scalar(pVec[6], pVec[8], 0, INFINITY));
    choice[idx[j]] = (dt1s2[j] < dt2s2[j]) ? 1 : 2;
    winDT[idx[j]]  = (dt1s2[j] < dt2s2[j]) ? (dt1s2[j] + pVec[9]) : (dt2s2[j] + pVec[9]);
  }
  
  arma::vec DT = winDT + pVec[7]; // Add t_delay; Note t_ND has yet been added 
  arma::mat out = arma::join_horiz(choice, DT);
  return out;
}

//' @rdname rplba
//' @export
// [[Rcpp::export]]
arma::mat rplba_omp(int n, arma::vec pVec) {
  /* A  b muv1  muv2  t_ND  muw1  muw2 t_delay  sv  swt; c(0, 1, 2, 3, 4)
     0  1    2     3     4     5     6       7   8    9; c(5, 6, 7) */
  arma::vec x1(n), x2(n), v1(n), v2(n), dt1(n), dt2(n);  
  arma::vec choice(n), winDT(n), undone(n);
  
#ifdef _OPENMP
#pragma omp parallel for default(shared) firstprivate(pVec)
#endif
  for(int i=0; i<n; i++) 
  { // Stage 1 LBA
    v1[i]  = rtn_scalar(pVec[2], pVec[8], 0, INFINITY) ; // acc1 (X)
    v2[i]  = rtn_scalar(pVec[3], pVec[8], 0, INFINITY) ; // acc2 (O)
    x1[i]  = Rf_runif(0, pVec[0]);
    x2[i]  = Rf_runif(0, pVec[0]);
    dt1[i] = (pVec[1] - x1[i]) / v1[i];
    dt2[i] = (pVec[1] - x2[i]) / v2[i];
    choice[i] = (dt1[i] < dt2[i]) ? 1 : 2;
    winDT[i]  = (dt1[i] < dt2[i]) ? dt1[i] : dt2[i];
    undone[i] = (winDT[i] <= pVec[9]) ? false : true;
  }
  
  int n2         = arma::accu(undone);
  arma::uvec idx = find(undone);
  arma::vec x1s2 = x1.elem(idx);  arma::vec x2s2 = x2.elem(idx);
  arma::vec v1s2 = v1.elem(idx);  arma::vec v2s2 = v2.elem(idx); 
  arma::vec dt1s2(n2), dt2s2(n2);
  
#ifdef _OPENMP
#pragma omp parallel for default(shared) firstprivate(pVec)
#endif
  for (int j=0; j<n2; j++) 
  { // Stage 2 LBA
    dt1s2[j] = (pVec[1] - (x1s2[j] + pVec[9]*v1s2[j])) / 
      (rtn_scalar(pVec[5], pVec[8], 0, INFINITY));
    dt2s2[j] = (pVec[1] - (x2s2[j] + pVec[9]*v2s2[j])) / 
      (rtn_scalar(pVec[6], pVec[8], 0, INFINITY));
    choice[idx[j]] = (dt1s2[j] < dt2s2[j]) ? 1 : 2;
    winDT[idx[j]]  = (dt1s2[j] < dt2s2[j]) ? (dt1s2[j] + pVec[9]) : (dt2s2[j] + pVec[9]);
  }
  
  arma::vec DT = winDT + pVec[7]; // Add t_delay; Note t_ND has yet been added 
  arma::mat out = arma::join_horiz(choice, DT);
  return out;
}

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


