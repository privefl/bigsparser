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

  ofstream myFile(filename, ios::out | ios::binary);

  const int    * pI = i.begin();
  const double * pX = x.begin();
  int K = i.size();

  for (int k = 0; k < K; k++, pI++, pX++) {
    myFile.write(reinterpret_cast<const char*>(pI), 4);
    myFile.write(reinterpret_cast<const char*>(pX), 8);
  }

  myFile.close();
}

/******************************************************************************/
