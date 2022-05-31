#ifndef SFBM_H
#define SFBM_H

/******************************************************************************/

#ifndef STRICT_R_HEADERS
#define STRICT_R_HEADERS
#endif

#include <mio/mmap.hpp>
#include <system_error> // for std::error_code
#include <Rcpp.h>

/******************************************************************************/

class SFBM {
public:
  SFBM(std::string path, int n, int m, const std::vector<size_t>& p,
       const std::vector<int>& first_i)
  : n(n), m(m), p(p), first_i(first_i) {

    std::error_code error;
    this->ro_mmap.map(path, error);
    if (error) Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());

    this->data = reinterpret_cast<const double*>(ro_mmap.data());

    this->compact = first_i.size() > 0;
  }

  const double* i_x() { return reinterpret_cast<const double*>(ro_mmap.data()); }
  size_t nrow() const { return n; }
  size_t ncol() const { return m; }
  bool is_compact() const { return compact; }

  template<class C>
  C prod(const C& x) {

    C res(n);
    for (int i = 0; i < n; i++) res[i] = 0;

    for (int j = 0; j < m; j++) {

      double x_j = x[j];

      if (x_j != 0) {

        if (compact) {

          size_t lo = p[j];
          size_t up = p[j + 1];
          int i = first_i[j];

          for (size_t k = lo; k < up; k++, i++)
            res[i] += data[k] * x_j;

        } else {

          size_t lo = 2 * p[j];
          size_t up = 2 * p[j + 1];

          for (size_t k = lo; k < up; k += 2) {
            int    ind = data[k];
            double val = data[k + 1];
            res[ind] += val * x_j;
          }

        }
      }
    }

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

    if (compact) {

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

    } else {

      size_t lo = 2 * p[j];
      size_t up = 2 * p[j + 1];

      size_t k = lo;

      if (up >= (lo + 8)) {
        for (; k <= (up - 8); k += 8) {  // unrolling optimization
          cp += (data[k + 1] * x[data[k]] + data[k + 3] * x[data[k + 2]]) +
            (data[k + 5] * x[data[k + 4]] + data[k + 7] * x[data[k + 6]]);
        }
      }
      for (; k < up; k += 2) {
        int    ind = data[k];
        double val = data[k + 1];
        cp += val * x[ind];
      }

    }

    return cp;
  }

  template<class C>
  C& incr_mult_col(int j, C& x, double coef) {

    if (compact) {

      size_t lo = p[j];
      size_t up = p[j + 1];
      int i = first_i[j];

      for (size_t k = lo; k < up; k++, i++)
        x[i] += data[k] * coef;

    } else {

      size_t lo = 2 * p[j];
      size_t up = 2 * p[j + 1];

      for (size_t k = lo; k < up; k += 2) {
        int    ind = data[k];
        double val = data[k + 1];
        x[ind] += val * coef;
      }

    }

    return x;
  }

protected:
  mio::mmap_source ro_mmap;
  const double * data;
  int n;
  int m;
  std::vector<size_t> p;
  std::vector<int> first_i;
  bool compact;
};

/******************************************************************************/

#endif // SFBM_H
