/******************************************************************************/

#include <bigsparser/SFBM.h>

/******************************************************************************/

// [[Rcpp::export]]
SEXP getXPtrSFBM(std::string path, int n, int m, std::vector<size_t> p) {

  if (int(p.size()) != (m + 1)) Rcpp::stop("This is a bug.");

  // http://gallery.rcpp.org/articles/intro-to-exceptions/
  try {
    // Create a pointer to an FBM object and wrap it as an external pointer
    Rcpp::XPtr<SFBM> ptr(new SFBM(path, n, m, p), true);
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
