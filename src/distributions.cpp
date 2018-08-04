#include <cpda.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <boost/math/distributions/inverse_gaussian.hpp> // for inverse_gaussian_distribution
using boost::math::inverse_gaussian; // typedef provides default type is double.
using boost::math::inverse_gaussian_distribution; // for inverse gaussian distribution.

using namespace Rcpp;

void set_seed(unsigned int seed) {
  Environment base_env("package:base");
  Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

/* These three functions are to draw random number, ie rexGaussian */ 
struct exG_params {
  double mu, sigma, nu, prob;
};

inline double pexGAUS_(double q, double mu, double sigma, double nu, 
  bool lt = true, bool lp = false) {
  
  double cdf, term1, term2, term3, term4, term5, z, sigma2, out;
  
  if (nu > 0.05 * sigma) {
    sigma2 = sigma*sigma;
    z = q - mu - (sigma2/ nu);
    term1 = R::pnorm5( (q - mu)/sigma, 0, 1, true, false);
    term2 = R::pnorm5(z/sigma, 0, 1, true, false);
    term3 = mu + sigma2/nu;
    term4 = mu*mu + 2.0*q*sigma2/nu;
    term5 = std::exp( (term3*term3 - term4) / (2.0*sigma2) );
    cdf = lt ? term1 - term2 * term5 : 1 - (term1 - term2 * term5);
  } else {
    cdf = lt ? R::pnorm5(q, mu, sigma, true, false) : (1 - R::pnorm5(q, mu, sigma, true, false));
  }
  
  out = !lp ? cdf : std::log(cdf);
  
  return out;
}

double pexGaussian_function (double q, void *params) {
  struct exG_params *p 
  = (struct exG_params *) params;
  
  double mu = p->mu;
  double sigma = p->sigma;
  double nu = p->nu;
  double prob = p->prob;
  double out = pexGAUS_(q, mu, sigma, nu, true, false) - prob;
  
  return out;
}

inline double rtnorm0(const double l, const double u) {
  // Accept-Reject Algorithm 0; Naive method A-R method
  bool invalid = true;
  double z;
  while (invalid) {
    z = R::rnorm(0.0, 1.0);
    if (z <= u && z >= l) break;
  }
  return z;
}

inline double rtnorm1(const double l, const double u) {
  // Algorithm 1; 'expl'; use when lower > mean; upper = INFINITY; p 122, right
  bool invalid = true;
  double z, r, num; // a stands for alphaStar in Robert (1995)
  double a = 0.5 * (std::sqrt(l*l + 4.0) + l);
  
  while (invalid) {
    z   = (-1.0/a)*std::log(R::runif(0.0, 1.0)) + l; // control lower boundary
    num = R::runif(0.0, 1.0);
    r   = std::exp(-0.5 * (z - a)*(z - a));
    if (num <= r && z <= u) break;
  }
  return z ;
}

inline double rtnorm2(const double l, const double u) {
  // Algorithm 2; 'expu'; use when upper < mean; lower = -INFINITY.
  bool invalid = true;
  double z, r, num;
  double a = 0.5 * (std::sqrt(u*u + 4.0) - u);
  
  while (invalid) {
    z   = (-1.0/a)*std::log(R::runif(0.0, 1.0)) - u; // control lower boundary
    num = R::runif(0.0, 1.0);
    r   = std::exp(-0.5 * (z - a)*(z - a));
    if (num <= r && z <= -l) break;
  }
  return -z;  // note the negative
}

inline double rtnorm3(const double l, const double u) {
  // Algorithm 3; page 123. 2.2. Two-sided truncated normal dist.
  bool invalid = true;
  double z, r, num; // a stands for alphaStar in Robert (1995)
  double l2 = l*l;
  double u2 = u*u;
  
  while (invalid) {
    z = R::runif(l, u) ;
    if (l > 0) {
      r = std::exp(0.5 * (l2 - z*z));
    } else if (u < 0) {
      r = std::exp(0.5 * (u2 - z*z));
    } else  {
      r = std::exp( -0.5 * z * z ) ;
    }
    num = R::runif(0.0, 1.0);
    if (num <= r) break;
  }
  return z ;
}

// [[Rcpp::export]]
double rtn_scalar(double mean, double sd, double l, double u) {
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
  return(z*sd + mean);
}

double dtn_scalar(double x, double mean, double sd, double lower, double upper, 
  bool lp) {
  
  double out, numer, denom;
  if ((x >= lower) && (x <= upper)) {
    denom = R::pnorm(upper, mean, sd, true, false) - 
            R::pnorm(lower, mean, sd, true, false);
    numer = R::dnorm(x, mean, sd, lp);
    out = lp ? (numer - std::log(denom)) : (numer/denom);
  } else {
    out = lp ? -INFINITY : 0;
  }
  return(out);
}

//' @export
// [[Rcpp::export]]
arma::vec dtn_(arma::vec x, double mean, double sd, double lower, double upper, 
  bool lp) {

  double numer, denom;
  arma::vec out(x.n_elem);
  
  for (size_t i = 0; i < x.n_elem; i++) {
    if ((x(i) >= lower) && (x(i) <= upper)) {
      denom = R::pnorm(upper, mean, sd, true, false) - 
        R::pnorm(lower, mean, sd, true, false);
      numer = R::dnorm(x(i), mean, sd, lp);
      out(i) = lp ? (numer - std::log(denom)) : (numer/denom);
    } else {
      out(i) = lp ? -INFINITY : 0;
    }
  }
  return out;
}

//' @export
// [[Rcpp::export]]
arma::vec dtn(arma::vec x, arma::vec mean, arma::vec sd, arma::vec lower,
  arma::vec upper, bool log = false) {
  unsigned int n = x.n_elem;
  if (mean.n_elem == 1) mean = arma::repmat(mean, n, 1);
  if (sd.n_elem == 1) sd = arma::repmat(sd, n, 1);
  if (lower.n_elem == 1) lower = arma::repmat(lower, n, 1);
  if (upper.n_elem == 1) upper = arma::repmat(upper, n, 1);
  arma::vec out(x.n_elem);
  
  for (size_t i = 0; i < n; i++)
    out(i) = dtn_scalar(x(i), mean(i), sd(i), lower(i), upper(i), log);
  return out;
}


double ptn_scalar(double q, double mean, double sd, double lower, double upper,
  bool lt, bool lp) {
  
  double out, numer, denom, qtmp;
  if (lt) {
    out = (q < lower) ? 0 : 1;
  } else {
    out = (q < lower) ? 1 : 0;
  }
  if ((q >= lower) && (q <= upper)) {
    // 4th arg: lower.tail (lt)=1; 5th arg: log.p (lg)=0
    denom = R::pnorm(upper, mean, sd, true, false) -
            R::pnorm(lower, mean, sd, true, false);
    qtmp  = lt ? (R::pnorm(q, mean, sd, true, false) - R::pnorm(lower, mean, sd, true, false)) :
                 (R::pnorm(upper, mean, sd, true, false) - R::pnorm(q, mean, sd, true, false));
    out  = lp ? (std::log(qtmp)-std::log(denom)) : (qtmp/denom);
  }
  return out;
}

//' @export
// [[Rcpp::export]]
arma::vec ptn_(arma::vec q, double mean, double sd, double lower, double upper,
  bool lt, bool lp) {
  
  arma::vec out(q.n_elem);
  double numer, denom, qtmp;
  
  for (size_t i = 0; i < q.n_elem; i++) {
    if (lt) {
      out(i) = (q(i) < lower) ? 0 : 1;
    } else {
      out(i) = (q(i) < lower) ? 1 : 0;
    }
    
    if ((q(i) >= lower) && (q(i) <= upper)) {
      denom = R::pnorm(upper, mean, sd, true, false) -
              R::pnorm(lower, mean, sd, true, false);
      qtmp  = lt ? (R::pnorm(q(i), mean, sd, true, false) - R::pnorm(lower, mean, sd, true, false)) :
                   (R::pnorm(upper, mean, sd, true, false) - R::pnorm(q(i), mean, sd, true, false));
      out(i)  = lp ? (std::log(qtmp)-std::log(denom)) : (qtmp/denom);
    }
  }
  
  return out;
}

//' @export
// [[Rcpp::export]]
arma::vec ptn(arma::vec q, arma::vec mean, arma::vec sd, arma::vec lower,
  arma::vec upper, bool lt = true, bool lp = false) {

  unsigned int n = q.n_elem;
  if (mean.n_elem == 1) mean = arma::repmat(mean, n, 1);
  if (sd.n_elem == 1) sd = arma::repmat(sd, n, 1);
  if (lower.n_elem == 1) lower = arma::repmat(lower, n, 1);
  if (upper.n_elem == 1) upper = arma::repmat(upper, n, 1);
  arma::vec out(n);
  
  for(size_t i = 0; i < n; i++)
    out(i) = ptn_scalar(q(i), mean(i), sd(i), lower(i), upper(i), lt, lp);
  return out;
}

//' @rdname dtnorm
//' @export
// [[Rcpp::export]]
arma::vec rtnorm(unsigned int n, double mean, double sd, double lower,
  double upper) {
  arma::vec out(n);
  for (size_t i = 0; i < n; i++)
    out(i) = rtn_scalar(mean, sd, lower, upper);
  return out;
}



//' @export
// [[Rcpp::export]]
arma::vec dexGAUS(arma::vec x, arma::vec mu, arma::vec sigma,
   arma::vec tau, bool log = false) {
  if (any(sigma < 0)) stop("sigma must be greater than 0 \n");
  if (any(tau < 0)) stop("tau must be greater than 0\n");

  double z, term1, term2;
  unsigned int n = x.n_elem;
  if (mu.n_elem == 1) mu = arma::repmat(mu, n, 1);
  if (sigma.n_elem == 1) sigma = arma::repmat(sigma, n, 1);
  if (tau.n_elem == 1) tau = arma::repmat(tau, n, 1);
  arma::vec out(n);

  for (size_t i = 0; i < n; i++) {
    if ( tau(i) > 0.05 * sigma(i) ) {
      z = x(i) - mu(i) - (sigma(i)*sigma(i) / tau(i));
      term1 = -std::log(tau(i)) - ( (z + (sigma(i)*sigma(i)) / (2*tau(i))) / tau(i) );
      term2 = R::pnorm5(z/sigma(i), 0, 1, true, true);
      out(i) = log ? term1 + term2 : std::exp(term1 + term2);
    } else {
      out(i) = R::dnorm4(x(i), mu(i), sigma(i), log);
    }
  }
  return out;
}

//' @export
// [[Rcpp::export]]
arma::vec pexGAUS(arma::vec q, arma::vec mu, arma::vec sigma,
  arma::vec tau, bool lt = true, bool lp = false) {
  if (any(sigma < 0)) stop("sigma must be greater than 0 \n");
  if (any(tau < 0)) stop("tau must be greater than 0\n");
  
  double cdf, term1, term2, term3, term4, term5, z, sigma2;
  
  unsigned int n = q.n_elem;
  if (mu.n_elem == 1) mu = arma::repmat(mu, n, 1);
  if (sigma.n_elem == 1) sigma = arma::repmat(sigma, n, 1);
  if (tau.n_elem == 1) tau = arma::repmat(tau, n, 1);
  arma::vec out(n);
  
  for (size_t i = 0; i < n; i++) {
    if ( tau(i) > 0.05 * sigma(i) ) {
      
      z = q(i) - mu(i) - (sigma(i)*sigma(i) / tau(i));
      sigma2 = sigma(i)*sigma(i);
      
      term1 = R::pnorm5( (q(i) - mu(i))/sigma(i), 0, 1, true, false);
      term2 = R::pnorm5(z/sigma(i), 0, 1, true, false);
      term3 = mu(i) + sigma2/tau(i);
      term4 = mu(i)*mu(i) + 2.0*q(i)*sigma2/tau(i);
      term5 = std::exp( (term3*term3 - term4) / (2.0*sigma2) );
      cdf = lt ? term1 - term2 * term5 : 1 - (term1 - term2 * term5);
    } else {
      cdf = lt ? R::pnorm5(q(i), mu(i), sigma(i), true, false) :
      1 - R::pnorm5(q(i), mu(i), sigma(i), true, false);
    }
    
    out(i) = !lp ? cdf : std::log(cdf);
  }
  
  return out;
}

//' @export
// [[Rcpp::export]]
arma::vec qexGAUS(arma::vec p, arma::vec mu, arma::vec sigma,
    arma::vec tau, bool lt = true, bool lp = false) {
    
    if (arma::any(sigma < 0)) stop("sigma must be greater than 0 \n");
    if (arma::any(tau < 0)) stop("tau must be greater than 0\n");
    if (lp == true) p = arma::exp(p); 
    if (lt == false) p = 1 - p;
    if (arma::any(p < 0) || arma::any(p > 1)) stop("p must be between 0 and 1\n");     
    
    unsigned int n = p.n_elem;
    double x_lo, x_hi, r;
    
    if (mu.n_elem == 1) mu = arma::repmat(mu, n, 1);
    if (sigma.n_elem == 1) sigma = arma::repmat(sigma, n, 1);
    if (tau.n_elem == 1) tau = arma::repmat(tau, n, 1);
    arma::vec out(n);
    
    for (size_t i = 0; i < n; i++) {
      bool sw1 = pexGAUS_(mu(i), mu(i), sigma(i), tau(i), true, false) < p(i);
      if (sw1) {
        x_lo = mu(i);
        x_hi = mu(i) + sigma(i);
        unsigned int j = 1;
        while (pexGAUS_(x_hi, mu(i), sigma(i), tau(i), true, false) < p(i)) {
          x_hi = mu(i) + (double)j*sigma(i);
          j++;
        } 
      } else {
        x_lo = mu(i) - sigma(i);
        x_hi = mu(i);
        unsigned int j = 1;
        while (pexGAUS_(x_lo, mu(i), sigma(i), tau(i), true, false) > p(i)) {
          x_lo = mu(i) - (double)j*sigma(i);
          j++;
        }
      }
      
      int status;
      unsigned int iter = 0, max_iter = 100;
      const gsl_root_fsolver_type *T;
      gsl_root_fsolver *s;
      gsl_function F;
      struct exG_params params = {mu(i), sigma(i), tau(i), p(i)};
      
      F.function = &pexGaussian_function;
      F.params = &params;
      T = gsl_root_fsolver_brent;
      s = gsl_root_fsolver_alloc (T);
      gsl_root_fsolver_set (s, &F, x_lo, x_hi);
      
      do {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
      } while (status == GSL_CONTINUE && iter < max_iter);
      
      gsl_root_fsolver_free (s);
      out(i) = r; 
    }

    return out;
}

//' @export
// [[Rcpp::export]]
arma::vec dinvGAUS(arma::vec x, arma::vec mu, arma::vec lambda, bool log = false) { 
  
  if (any(mu <= 0)) stop("mu must be greater than 0 \n");
  if (any(lambda <= 0)) stop("lambda must be greater than 0\n");

  double den;
  unsigned int n = x.n_elem;
  if (mu.n_elem == 1) mu = arma::repmat(mu, n, 1);
  if (lambda.n_elem == 1) lambda = arma::repmat(lambda, n, 1);
  arma::vec out(n);
  
  for (size_t i = 0; i < n; i++) {
    inverse_gaussian w11(mu(i), lambda(i));
    den = pdf(w11, x(i));
    out(i) = log ? std::log(den) : den;
  }
  
  return out;
}

//' @export
// [[Rcpp::export]]
arma::vec pinvGAUS(arma::vec q, arma::vec mu, arma::vec lambda, bool lt = true, 
  bool lp = false) {
  
  if (any(mu <= 0)) stop("mu must be greater than 0 \n");
  if (any(lambda <= 0)) stop("lambda must be greater than 0\n");

  double den, tmp;
  unsigned int n = q.n_elem;
  if (mu.n_elem == 1) mu = arma::repmat(mu, n, 1);
  if (lambda.n_elem == 1) lambda = arma::repmat(lambda, n, 1);
  arma::vec out(n);
  
  for (size_t i = 0; i < n; i++) {
    inverse_gaussian w11(mu(i), lambda(i));
    tmp = cdf(w11, q(i));
    den = lt ? tmp : 1 - tmp;   
    out(i) = !lp ? den : std::log(den);
  }
  
  return out;
}

//' @export
// [[Rcpp::export]]
arma::vec qinvGAUS(arma::vec p, arma::vec mu, arma::vec lambda, bool lt = true, 
  bool lp = false) {
  
  if (any(mu <= 0)) stop("mu must be greater than 0 \n");
  if (any(lambda <= 0)) stop("lambda must be greater than 0\n");
  if (lp == true) p = arma::exp(p); 
  if (lt == false) p = 1 - p;
  if (arma::any(p < 0) || arma::any(p > 1)) stop("p must be between 0 and 1\n");     
  
  double den, tmp;
  unsigned int n = p.n_elem;
  if (mu.n_elem == 1) mu = arma::repmat(mu, n, 1);
  if (lambda.n_elem == 1) lambda = arma::repmat(lambda, n, 1);
  arma::vec out(n);
  
  for (size_t i = 0; i < n; i++) {
    inverse_gaussian w11(mu(i), lambda(i));
    out(i) = quantile(w11, p(i));
  }
  
  return out;
}


//' @export
// [[Rcpp::export]]
arma::vec rinvGAUS_(unsigned int n, double mu, double lambda) {
  if (mu <= 0) stop("mu must be greater than 0 \n");
  if (lambda <= 0) stop("lambda must be greater than 0\n");
  double a, c, v, x, half_mu_lambda, mu2;

  half_mu_lambda = 0.5*mu/lambda; // C in Michael, Schucany & Hass's FORTRAN subroutine
  mu2 = mu*mu;
  a = mu*half_mu_lambda;
  c = 4.0*mu*lambda;
  arma::vec out = arma::randn(n);
  
  for (size_t i = 0; i < n; i++) {
      v = out(i) * out(i);	// Chi-square with 1 df
      x = mu + a*v - half_mu_lambda*sqrt(c*v+mu2*v*v);	// Smallest root
      out(i) = (R::runif(0, 1) < (mu/(mu+x))) ? x : mu2/x;	// Pick x with prob mu/(mu+x), else d/x;
      if (out(i) < 0.0) v = x;
  }
  return out;
}

//' @export
// [[Rcpp::export]]
arma::vec rinvGAUS(unsigned int n, arma::vec mu, arma::vec lambda) {
  if (any(mu <= 0)) stop("mu must be greater than 0 \n");
  if (any(lambda <= 0)) stop("lambda must be greater than 0\n");

  double a, c, v, x, half_mu_lambda, mu2;
  arma::vec out = arma::randn(n);

  for (size_t i = 0; i < n; i++) {
    half_mu_lambda = 0.5*mu(i)/lambda(i); 
    mu2 = mu(i)*mu(i);
    a = mu(i)*half_mu_lambda;
    c = 4.0*mu(i)*lambda(i);
    
    v = out(i) * out(i);	// Chi-square with 1 df
    x = mu(i) + a*v - half_mu_lambda*sqrt(c*v+mu2*v*v);	// Smallest root
    out(i) = (R::runif(0, 1) < (mu(i)/(mu(i)+x))) ? x : mu2/x;	// Pick x with prob mu/(mu+x), else d/x;
  }

  return out;
}
