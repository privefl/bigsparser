#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sum_overflow(IntegerVector x) {
  return Rcpp::sum(x);
}

/*** R
sum_overflow(seq_len(1e6))  # 1784293664 = (500000500000 %% 2^32)
sum(seq_len(1e6))  # 500000500000
*/
