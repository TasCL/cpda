#include <cpda.hpp>
#include <random>

// Protect against compilers without OpenMP; e.g., OS X  clang
#ifdef _OPENMP
#include <omp.h>
#endif

arma::vec pmax(arma::vec v, double min) {
  for (arma::vec::iterator it=v.begin(); it!=v.end(); it++)
  {
    if (*it < min) *it = min ;
  }
  return v ;
}
arma::vec getVec(double *x, int *nx) {
  arma::vec out(*nx);
  for(int i=0; i<*nx; i++) { out[i]=*(x+i); }
  return out;
}
double gaussian(double y, arma::vec yhat, double h) {
  // standard gaussian kernel mean=0; sigma==1
  double x;
  int ns = yhat.n_elem;
  arma::vec result(ns);
  
  for(arma::vec::iterator it=yhat.begin(); it!=yhat.end(); ++it)
  {
    int i = std::distance(yhat.begin(), it);
    x = (y - *it)/h;  // z / h
    // (1/h) * K(z/h); K_h(z)
    result[i] = ( (1/(sqrt(2*arma::datum::pi))) * exp( -pow(x,2) / 2 ) ) / h;
  }
  // (1/N_s) * sigma K_h (x-x.tidle_j)
  return ( arma::accu(result) / ns);
}
double gaussian_omp(double y, arma::vec yhat, double h) {
  // standard gaussian kernel mean=0; sigma==1
  double x;
  int ns = yhat.n_elem;
  arma::vec result(ns);
  
#ifdef _OPENMP
#pragma omp parallel for default(shared) firstprivate(y, yhat, h)  
#endif
  for(int j=0; j < ns; j++)
  {
    x = (y - yhat[j])/h;  
    result[j] = ( (1/(sqrt(2*arma::datum::pi))) * exp( -pow(x,2) / 2 ) ) / h;
  }
  
  return ( arma::accu(result) / ns);
}

//' A simple and fast quantile calculator
//'
//' A C++ quantile function.
//'
//' @param y a data vector
//' @param q nth quantile. Enter proportion, such as .25 or .75.
//' @examples
//' y <- rnorm(100)
//' q <- cquantile(y, .25) 
//' @export
// [[Rcpp::export]]
double cquantile(arma::vec y, double q) {
  arma::vec sy = sort(y) ;
  int ny = sy.n_elem ;
  int nth = ny * (q - 1e-9) ;
  return sy(nth);
}

//' Silverman's Rule of Thumb Bandwidth for Kernel Density Estimation
//'
//' A C++ version of Silverman's rule-of-thumb bandwidth. This is similar
//' with R's \code{bw.nrd0(x)}
//'
//' @param y a data vector
//' @param m a multiplier to adjust the SROT proportionally.
//' 
//' @seealso
//' \code{\link{bw.nrd}}, \code{\link{bandwidth.nrd}}.
//' @export
//' @examples
//' data(lba)
//' h <-cpda::bwNRD0(plba$DT1, 0.8)
//' 
// [[Rcpp::export]]
double bwNRD0(arma::vec y, double m) {
  int n = y.n_elem ;
  double out ;
  if(n < 2) {
    // fprintf(stderr, "need at least 2 data points\n");
    out = 0.9 * y[0] * pow(n, -0.2) ;
  } else {
    double q25 = cquantile(y, 0.25) ;
    double q75 = cquantile(y, 0.75) ;
    // double h   = (q75-q25)/1.34 ;
    double h   = (q75-q25);
    double sd  = arma::stddev(y) ;
    double IS  = std::min(h, sd);
    out = 0.9 * IS * pow(n, -0.2) ;
  }
  return m*out ;
}

struct genInt {
  int x ;
  genInt() {x=0;}
  int operator()() {return x++;}
} ; // class generator:

std::vector<int> generateIntVec (int n) {
  // generate an integer vector based on C index from 0 to n-1 
  std::vector<int> out(n); // a empty vector with n integer mem space.
  genInt generator; // 
  std::generate(out.begin(), out.end(), generator);
  return out;
}

std::vector<int> shuffle(int n) {
  std::vector<int> out = generateIntVec(n); 
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle(out.begin(), out.end(), std::default_random_engine(seed));
  return out ;
} 
