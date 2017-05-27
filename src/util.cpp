#include <cpda.hpp>
#include <random>

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
  int n = y.n_elem;
  return m*0.9*std::min((cquantile(y, .75) - cquantile(y, .25)),
  arma::stddev(y)) * std::pow((double)n, -.2);
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

//' @export
// [[Rcpp::export]]
arma::uvec histc(arma::vec x, arma::vec edge) {
  arma::uvec hc = arma::histc(x, edge) ;
  // arma::vec bincount  = arma::conv_to<arma::vec>::from(hc);
  return hc ;
}

//' @export
// [[Rcpp::export]]
arma::vec histd(arma::vec yhat, arma::vec z, int n) {
  double dt, total_count, unit;
  arma::vec bin_edges, bc, pdf, out;
  // 1. Use the range of empirical data to define min and max
  // 2. Get a linear spaced 2^p numbers within [min-3*h, max+3*h] 
  dt = z[1] - z[0];
  bc = arma::conv_to<arma::vec>::from(arma::histc(yhat, getEdges(z))); // 1025
  // If this is a defective PDF, the user must tell 'histd' the total count.
  if(n==0) {total_count = arma::accu(bc);} else {total_count=n;}
  unit = dt * total_count;
  arma::vec unitVec(bc.n_elem); unitVec.fill(unit);
  pdf = bc / unitVec;
  out = pdf.rows(0, pdf.size()-2); // 0 to 1023
  return out;
}

//' Calculate Histogram Edges 
//'
//' This is an internal function. The user may not use it.
//' 
//' @param z a grid a scalar, usually created by 
//' \code{z = arma::linspace<arma::vec>(z0, z1, 1<<(int)p);}
//' 
//' @examples
//' set.seed(123)
//' dat1 <- stats::rnorm(1e2)
//' h  <- bw.nrd0(dat1)
//' z0 <- min(dat1) - 3*h 
//' z1 <- max(dat1) + 3*h
//' ngrid <- 2^10
//' dt <- (z1 - z0) / (ngrid - 1)
//' z  <- z0 + (0:(ngrid-1)) * dt
//' ## Same as using seq function
//' ## seq(z0, z1, length.out=ngrid)
//' 
//' binedge1 <- c(z - 0.5*dt, z[ngrid] + 0.5*dt)
//' binedge2 <- as.vector(cpda::getEdges(z))
//' all.equal(binedge1, binedge2)
//' mean(binedge - as.vector(tmp))
//' ## [1] 2.640334e-17
//' 
//' ## 4.00 vs. 14.1 microseconds
//' ## library(microbenchmark)
//' ## res <- microbenchmark(
//' ##    binedge1 <- c(z - 0.5*dt, z[ngrid] + 0.5*dt),
//' ##    binedge2 <- as.vector(cpda::getEdges(z)),
//' ##   times=10L)
//' @export
// [[Rcpp::export]]
arma::vec getEdges(arma::vec z) {
  double halfdt = (z[1] - z[0])/2;   
  arma::vec dtVec(z.n_elem);
  dtVec.fill(halfdt);
  arma::vec term1 = z - dtVec;
  arma::vec term2(1);
  double last     = z[z.n_elem - 1] + halfdt;
  term2.fill(last);
  return arma::join_cols(term1, term2) ;
}

//' @export
// [[Rcpp::export]]
arma::vec getFilter(double m, double M, double h, double p) {
  // See gaussian filter equation in https://en.wikipedia.org/wiki/Gaussian_filter
  double N_grid  = std::pow(2.0, p);
  double tmp0    = arma::datum::pi * N_grid/(M-m) ;
  arma::vec tmp1 = arma::linspace<arma::vec>(0, 1, 1 + N_grid/2) ;
  
  arma::vec tmp0Vec(tmp1.n_elem);
  tmp0Vec.fill(tmp0);
  arma::vec freq = tmp0Vec % tmp1 ;
  arma::vec freq2= arma::pow(freq, 2.0) ; // s^2 on p17
  
  double h2      = std::pow(h, 2.0) ;
  arma::vec fil0 = arma::exp(-0.5 * h2 * freq2) ;
  arma::vec fil1 = arma::flipud(fil0.rows(1, (fil0.size() - 2)));
  
  arma::vec out  = arma::join_cols(fil0, fil1) ;
  return out ;
}

