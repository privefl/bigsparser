/******************************************************************************/

#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;
using namespace std;

/******************************************************************************/

// [[Rcpp::export]]
void write_indval(const char * filename,
                  const IntegerVector& i,
                  const NumericVector& x) {

  ofstream myFile(filename, ios::out | ios::app | ios::binary);

  size_t K = i.size();

  double dbl_i = 0, dbl_x = 0;
  const char * ptr_char_i = reinterpret_cast<const char*>(&dbl_i);
  const char * ptr_char_x = reinterpret_cast<const char*>(&dbl_x);

  for (size_t k = 0; k < K; k++) {
    dbl_i = i[k];
    myFile.write(ptr_char_i, 8);
    dbl_x = x[k];
    myFile.write(ptr_char_x, 8);
  }

  myFile.close();
}

/******************************************************************************/
