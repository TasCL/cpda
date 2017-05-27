#include <cpda.hpp>

void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

#define SQRT_2PI 2.5066282746310007e+0  /* sqrt(2 x pi) */
#define M_E	     2.7182818284590452354  /* e */

inline double rtnorm0(double l, double u) {
  // Accept-Reject Algorithm 0; Naive method A-R method
  int valid = 0;
  double z;
  while (valid == 0) {
    z = R::rnorm(0.0, 1.0);
    if (z <= u && z >= l) { break;}
  }
  return z ;
}

inline double rtnorm1(double l, double u) {
  // Algorithm 1; 'expl'; use when lower > mean; upper = INFINITY; p 122, right
  int valid = 0;
  double z, r, num; // a stands for alphaStar in Robert (1995)
  double a = 0.5 * (sqrt(l*l + 4.0) + l);
  
  while (valid == 0) {
    z   = (-1/a)*std::log(R::runif(0.0, 1.0)) + l; // control lower boundary
    num = R::runif(0.0, 1.0);
    r   = exp(-0.5 * (z - a)*(z - a));
    if (num <= r && z <= u) { break; }
  }
  return z ;
}

inline double rtnorm2(double l, double u) {
  // Algorithm 2; 'expu'; use when upper < mean; lower = -INFINITY.
  int valid = 0;
  double z, r, num;
  double a = 0.5 * (sqrt(u*u + 4.0) - u);
  
  while (valid == 0) {
    z   = (-1/a)*std::log(R::runif(0.0, 1.0)) - u; // control lower boundary
    num = R::runif(0.0, 1.0);
    r   = exp(-0.5 * (z - a)*(z - a));
    if (num <= r && z <= -l) { break; }
  }
  return -z;  // note the negative
}

inline double rtnorm3(double l, double u) {
  // Algorithm 3; page 123. 2.2. Two-sided truncated normal dist.
  int valid = 0;
  double z, r, num; // a stands for alphaStar in Robert (1995)
  double l2 = l*l;
  double u2 = u*u;
  
  while (valid == 0) {
    z = R::runif(l, u) ;
    if (l > 0) {
      r = std::exp(0.5 * (l2 - z*z));
    } else if (u < 0) {
      r = std::exp(0.5 * (u2 - z*z));
    } else  {
      r = std::exp( -0.5 * z * z ) ;
    }
    num = R::runif(0, 1);
    if (num <= r) { break; }
  }
  return z ;
}

inline double rtn_scalar(double mean, double sd, double l, double u) {
  double z, stdl, stdl2, stdu, stdu2, eq_a1, eq_a2; // Standardised lower and upper
  bool a0, a1, a2, a3;
  stdl  = (l - mean) / sd; // l == stdlower, u == stdupper
  stdu  = (u - mean) / sd;
  stdl2 = stdl*stdl;
  stdu2 = stdu*stdu;
  
  // Accept-Reject Algorithm 0;
  // Algorithm (1): Use Proposition 2.3 with only lower truncation. upper==INFINITY
  // rejection sampling with exponential proposal. Use if lower > mean
  // Algorithm (2): Use -x ~ N_+ (-mu, -mu^+, sigma^2) on page 123. lower==-INFINITY
  // rejection sampling with exponential proposal. Use if upper < mean.
  // Algorithm (3, else): rejection sampling with uniform proposal.
  // Use if bounds are narrow and central.
  a0 = (stdl < 0 && u == INFINITY) || (stdl == -INFINITY && stdu > 0) ||
    (std::isfinite(stdl) && std::isfinite(u) && stdl < 0 && stdu > 0 && (stdu - stdl) > SQRT_2PI);
  eq_a1 = stdl + (2.0 * std::sqrt(M_E) / (stdl + std::sqrt(stdl2 + 4.0))) *
    (std::exp( 0.25 * (2.0*stdl - stdl*std::sqrt(stdl2 + 4.0))));
  a1 = (stdl >= 0) && (stdu > eq_a1);
  eq_a2 = -stdu + (2.0 * std::sqrt(M_E) / (-stdu + std::sqrt(stdu2 + 4.0))) *
    (std::exp(0.25 * (2.0*stdu + stdu*std::sqrt(stdu2 + 4.0))));
  a2 = (stdu <= 0) && (-stdl > eq_a2);
  
  if (a0) {
    z = rtnorm0(stdl, stdu);
  } else if (a1) {
    z = rtnorm1(stdl, stdu);
  } else if (a2) {
    z = rtnorm2(stdl, stdu);
  } else {
    z = rtnorm3(stdl, stdu);
  }
  return z*sd + mean;
}

inline double dtn_scalar(double x, double mean, double sd, double lower, 
  double upper, int lp)
{
  double out, numer, denom;
  if ((x >= lower) && (x <= upper)) {
    // 4th arg: lower.tail (lt)=1; 5th arg: log.p (lg)=0
    denom = R::pnorm(upper, mean, sd, 1, 0) - R::pnorm(lower, mean, sd, 1, 0);
    numer = R::dnorm(x, mean, sd, lp);
    out = lp ? (numer - std::log(denom)) : (numer/denom);
  } else {
    out = lp ? -INFINITY : 0;
  }
  return(out);
}

double ptn_scalar(double q, double mean, double sd, double lower, double upper,
  int lt, int lp) {
  double out, numer, denom, qtmp;
  if (lt) {
    out = (q < lower) ? 0 : 1;
  } else {
    out = (q < lower) ? 1 : 0;
  }
  if ((q >= lower) && (q <= upper)) {
    denom = R::pnorm(upper, mean, sd, 1, 0) - R::pnorm(lower, mean, sd, 1, 0);
    qtmp  = lt ? (R::pnorm(q, mean, sd, 1, 0) - R::pnorm(lower, mean, sd, 1, 0)) :
      (R::pnorm(upper, mean, sd, 1, 0) - R::pnorm(q, mean, sd, 1, 0));
    out  = lp ? (std::log(qtmp)-std::log(denom)) : (qtmp/denom);
  }
  return(out);
}

//' Truncated Normal Distribution
//'
//' Three functions implementing random number generator, density function and
//' cumulative distribution function (\code{rtn}, \code{dtn}, and \code{ptn})
//' for truncation normal distribution.
//'
//' @param x,n,q x in \code{dtn} is a vector of quantiles; n in \code{rtn}
//' is number of observations. n must be a scalar. q is a vector of quantiles.
//' @param mean mean (must be scalar).
//' @param sd standard deviation (must be scalar).
//' @param lower lower truncation value (must be scalar).
//' @param upper upper truncation value (must be scalar).
//' @param lt lower.tail. a boolean switch; if TRUE (default) probabilities are
//' \code{P[X <= x]}, otherwise, \code{P[X > x]}.
//' @param lp log.p. a boolean switch; if TRUE (default is FALSE) probabilities p
//' are given as \code{log(p)}.
//' @return a column vector.
//' @examples
//' ## rtn example
//' dat1 <- rtn(1e5, 0, 1, 0, Inf)
//' ## dat2 <- msm::rtnorm(n, mean, sd, lower, upper)
//' ## den2 <- density(dat2)
//' hist(dat1, breaks="fd", freq=F)
//' ## lines(den2$x, den2$y, lwd=2.5)
//' ## res <- microbenchmark(
//' ##     rtn(n, mean, sd, lower, upper),
//' ##     msm::rtnorm(n, mean, sd, lower, upper))
//' ## print(res)
//'
//' ## dtn example
//' x <- seq(-5, 5, length.out=1e3)
//' dat1 <- dtn(x, 0, 1, -2, 2, 0)
//' ## dat2 <- msm::dtnorm(x, mean, sd, lower, upper, 0)
//' plot(x, dat1, type="l", lwd=2)
//' ## lines(x, dat2, lwd=2, lty="dashed", col="red")
//'
//' ## res <- microbenchmark(
//' ##     dtn(x, mean, sd, lower, upper, 0),
//' ##     msm::dtnorm(x, mean, sd, lower, upper, 0))
//' ## print(res)
//'
//' ## ptn example
//' x <- seq(-50, 10, length.out=1e3)
//' mean <- 0
//' sd <- 1
//' lower <- 0
//' upper <- 5
//' dat1 <- ptn(x, 0, 1, 0, 5, lp=TRUE)
//' ## dat2 <- msm::ptnorm(x, mean, sd, lower, upper, log.p=TRUE)
//' ## all(dat1[,1] == dat2)
//'
//' plot(x, log(dat1[,1]))
//' ## lines(x, log(dat2), col="red", lwd=2)
//' ## mtext("pnorm(x, log=TRUE)", adj = 0)
//' ## mtext("log(pnorm(x))", col = "red", adj = 1)
//' @export
// [[Rcpp::export]]
arma::vec dtn(arma::vec x, double mean, double sd, double lower, double upper,
  int lp=false) {
  if (upper < lower) {Rcpp::stop("'upper' must be greater than 'lower'.");}
  if (sd < 0)        {Rcpp::stop("'sd' must be greater than 0.\n");}
  if (sd==R_NegInf   || sd==R_PosInf)   {Rcpp::stop("'sd' must have a finite value.\n");}
  if (mean==R_NegInf || mean==R_PosInf) {Rcpp::stop("'mean' must have a finite value.\n");}
  
  int idx;
  arma::vec out(x.n_elem);
  
  for(arma::vec::iterator i = x.begin(); i < x.end(); i++) {
    idx = std::distance(x.begin(), i);
    out[idx] = dtn_scalar(*i, mean, sd, lower, upper, lp);
  }
  return out;
}

//' @rdname dtn
//' @export
// [[Rcpp::export]]
arma::vec rtn(int n, arma::vec mean, arma::vec sd, arma::vec lower,
  arma::vec upper) {
  arma::vec out(n);
  for(arma::vec::iterator i = out.begin(); i < out.end(); i++) {
    int idx = std::distance(out.begin(), i);
    *i = rtn_scalar(mean[idx], sd[idx], lower[idx], upper[idx]);
  }
  return out;
}

//' @rdname dtn
//' @export
// [[Rcpp::export]]
arma::vec ptn(arma::vec q, double mean, double sd, double lower, double upper,
  int lt=true, int lp=false) {
  if (upper < lower) {Rcpp::stop("'upper' must be greater than 'lower'.");}
  if (sd < 0)        {Rcpp::stop("'sd' must be greater than 0.\n");}
  if (sd==R_NegInf   || sd==R_PosInf)   {Rcpp::stop("'sd' must have a finite value.\n");}
  if (mean==R_NegInf || mean==R_PosInf) {Rcpp::stop("'mean' must have a finite value.\n");}
  
  int idx;
  arma::vec out(q.n_elem);
  
  for(arma::vec::iterator i = q.begin(); i < q.end(); i++) {
    idx = std::distance(q.begin(), i);
    out[idx] = ptn_scalar(*i, mean, sd, lower, upper, lt, lp);
  }
  return out;
  
}

//' @export
// [[Rcpp::export]]
arma::mat rlba(int n, double b, double A, arma::vec mean_v, arma::vec sd_v,
  double t0) {
  
  arma::vec RT(n), R(n);
  arma::vec finite_RT, finite_R;
  arma::uvec finite_idx;
  arma::mat out;
  double dt0, dt1;

  // inline double rtn_scalar(const double mean,  const double sd, const double l,
  //   const double u)
    
  for(int i=0; i<n; i++) {
    dt0 = (b - A*R::runif(0.0, 1.0)) / rtn_scalar(mean_v[0], sd_v[0], 0, INFINITY);
    dt1 = (b - A*R::runif(0.0, 1.0)) / rtn_scalar(mean_v[1], sd_v[1], 0, INFINITY);
    RT[i] = (dt0 < dt1) ? (dt0 + t0) : (dt1 + t0);
    R[i]  = (dt0 < dt1) ? 1 : 2;
  }
  
  if (RT.has_inf()) {
    Rcpp::Rcout << "Find Inf and remove them\n";
    finite_idx = arma::find_finite(RT); // trims off infinite values when parameters are
    finite_RT  = RT.elem(finite_idx);   // very improbable when fitting Bayesian models
    finite_R   = R.elem(finite_idx);
    out = arma::join_horiz(finite_RT, finite_R);
  } else {
    out = arma::join_horiz(RT, R);
  }
  return out;
}

//' @export
// [[Rcpp::export]]
arma::vec rlba_n1(int n, double b, double A, arma::vec mean_v, arma::vec sd_v,
  double t0) {

  arma::vec RT(n), R(n), finite_RT, finite_R, sRT0, out;
  arma::uvec finite_idx;
  double rt0, rt1;

  for(int i=0; i<n; i++) {
    rt0 = t0 + ((b - A*R::runif(0.0, 1.0)) / rtn_scalar(mean_v[0], sd_v[0], 0, INFINITY));
    rt1 = t0 + ((b - A*R::runif(0.0, 1.0)) / rtn_scalar(mean_v[1], sd_v[1], 0, INFINITY));
    RT[i] = (rt0 < rt1) ? rt0 : rt1;
    R[i]  = (rt0 < rt1) ? 1 : 2;
  }

  if (RT.has_inf()) {
    Rcpp::Rcout << "Find Inf and remove them\n";
    finite_idx = arma::find_finite(RT); // trims off infinite values when parameters are
    finite_RT  = RT.elem(finite_idx);   // very improbable when fitting Bayesian models
    finite_R   = R.elem(finite_idx);
    out  = finite_RT.elem(arma::find(finite_R==1));
  } else {
    out = RT.elem(arma::find(R==1));
  }
  return out;
}

 
