#include <cpda.hpp>

#include <boost/math/distributions/inverse_gaussian.hpp> // for inverse_gaussian_distribution
using boost::math::inverse_gaussian; // typedef provides default type is double.
using boost::math::inverse_gaussian_distribution; // for inverse gaussian distribution.

#include <boost/math/distributions/normal.hpp> // for normal_distribution
using boost::math::normal; // typedef provides default type is double.

#include <boost/array.hpp>
using boost::array;

#include <iostream>
using std::cout; using std::endl; using std::left; using std::showpoint; using std::noshowpoint;
#include <iomanip>
using std::setw; using std::setprecision;
#include <limits>
using std::numeric_limits;
#include <sstream>
using std::string;
#include <string>
using std::stringstream;

using namespace Rcpp;


//' @export
// [[Rcpp::export]]
void InvG(NumericVector x) {
  
  cout << "Example: Inverse Gaussian Distribution."<< endl;
  double tolfeweps = numeric_limits<double>::epsilon();

  int precision = 17; // traditional tables are only computed to much lower precision.
  cout.precision(17); // std::numeric_limits<double>::max_digits10; for 64-bit doubles.
  
  // Traditional tables and values.
  double step = 0.2; // in z
  double range = 4; // min and max z = -range to +range.
  
  // Construct a (standard) inverse gaussian distribution s
  inverse_gaussian w11(1, 1);
  // (default mean = units, and standard deviation = unity)
  cout << "(Standard) Inverse Gaussian distribution, mean = "<< w11.mean()
       << ", scale = " << w11.scale() << endl;
  
  cout << "Probability distribution function (pdf) values" << endl;
  cout << "  z " "      pdf " << endl;
  cout.precision(5);
  for (double z = (numeric_limits<double>::min)(); z < range + step; z += step)
  {
    cout << left << setprecision(3) << setw(6) << z << " "
         << setprecision(precision) << setw(12) << pdf(w11, z) << endl;
  }
  cout.precision(6); // default
  
}


