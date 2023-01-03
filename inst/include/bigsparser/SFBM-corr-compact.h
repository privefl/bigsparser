#ifndef SFBM_CORR_COMPACT_H
#define SFBM_CORR_COMPACT_H

/******************************************************************************/

#ifndef STRICT_R_HEADERS
#define STRICT_R_HEADERS
#endif

#include <mio/mmap.hpp>
#include <system_error> // for std::error_code
#include <Rcpp.h>

/******************************************************************************/

class SFBM_corr_compact {
public:
  SFBM_corr_compact(std::string path, int n, int m, const std::vector<size_t>& p,
                    const std::vector<int>& first_i)
  : n(n), m(m), p(p), first_i(first_i) {

    std::error_code error;
    this->ro_mmap.map(path, error);
    if (error)
      Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());

    this->data = reinterpret_cast<const int16_t*>(ro_mmap.data());
  }

  const int16_t* i_x() { return reinterpret_cast<const int16_t*>(ro_mmap.data()); }
  size_t nrow() const { return n; }
  size_t ncol() const { return m; }

  template<class C>
  C prod(const C& x) {

    C res(n);
    for (int i = 0; i < n; i++) res[i] = 0;

    for (int j = 0; j < m; j++) {

      double x_j = x[j];

      if (x_j != 0) {

        size_t lo = p[j];
        size_t up = p[j + 1];
        int i = first_i[j];

        for (size_t k = lo; k < up; k++, i++)
          res[i] += data[k] * x_j;

      }
    }

    for (int i = 0; i < n; i++) res[i] /= 32767;

    return res;
  }

  template<class C>
  C cprod(const C& x) {

    C res(m);

    for (int j = 0; j < m; j++)
      res[j] = dot_col(j, x);

    return res;
  }

  template<class C>
  double dot_col(int j, const C& x) {

    double cp = 0;

    size_t lo = p[j];
    size_t up = p[j + 1];

    size_t k = lo;
    int i = first_i[j];

    if (up >= (lo + 4)) {
      for (; k <= (up - 4); k += 4, i += 4) {  // unrolling optimization
        cp += (data[k] * x[i] + data[k + 1] * x[i + 1]) +
          (data[k + 2] * x[i + 2] + data[k + 3] * x[i + 3]);
      }
    }
    for (; k < up; k++, i++) cp += data[k] * x[i];

    return cp / 32767;
  }

  template<class C>
  C& incr_mult_col(int j, C& x, double coef) {

    double coef_scaled = coef / 32767;

    size_t lo = p[j];
    size_t up = p[j + 1];
    int i = first_i[j];

    for (size_t k = lo; k < up; k++, i++)
      x[i] += data[k] * coef_scaled;

    return x;
  }

protected:
  mio::mmap_source ro_mmap;
  const int16_t * data;
  int n;
  int m;
  std::vector<size_t> p;
  std::vector<int> first_i;
};

/******************************************************************************/

#endif // SFBM_CORR_COMPACT_H
