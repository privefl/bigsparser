/******************************************************************************/

#include <bigsparser/SFBM.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

/******************************************************************************/

SEXP getXPtrSFBM0(std::string path, int n, int m, const std::vector<size_t>& p,
                 const std::vector<int>& first_i) {

  if (int(p.size()) != (m + 1)) Rcpp::stop("This is a bug.");

  // http://gallery.rcpp.org/articles/intro-to-exceptions/
  try {
    // Create a pointer to an FBM object and wrap it as an external pointer
    Rcpp::XPtr<SFBM> ptr(new SFBM(path, n, m, p, first_i), true);
    // Return the external pointer to the R side
    return ptr;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("C++ exception (unknown reason)");
  }

  return R_NilValue;
}

/******************************************************************************/

// [[Rcpp::export]]
SEXP getXPtrSFBM(std::string path, int n, int m, const std::vector<size_t>& p) {

  std::vector<int> first_i(0);
  return getXPtrSFBM0(path, n, m, p, first_i);
}

// [[Rcpp::export]]
SEXP getXPtrSFBM_compact(std::string path, int n, int m,
                         const std::vector<size_t>& p,
                         const std::vector<int>& first_i) {

  return getXPtrSFBM0(path, n, m, p, first_i);
}

/******************************************************************************/
