/******************************************************************************/

#include <rmio/create-file.hpp>
#include <bigsparser/SFBM.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
void write_indval(std::string filename,
                  const IntegerVector& i,
                  const NumericVector& x,
                  size_t offset) {

  size_t K = i.size();

  if (offset == 0) {
    create_file(filename.c_str(), K, 16);
  } else {
    append_file(filename.c_str(), K, 16);
  }

  mio::mmap_sink rw_mmap;
  std::error_code error;
  rw_mmap.map(filename, error);
  if (error) Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());

  double * data = reinterpret_cast<double*>(rw_mmap.data());

  size_t k2 = 2 * offset;
  for (size_t k = 0; k < K; k++) {
    data[k2++] = i[k];
    data[k2++] = x[k];
  }
}

/******************************************************************************/
