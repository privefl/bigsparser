#ifndef SFBM_H
#define SFBM_H

/******************************************************************************/

#ifndef STRICT_R_HEADERS
#define STRICT_R_HEADERS
#endif

#include <mio/mmap.hpp>
#include <system_error> // for std::error_code
#include <Rcpp.h>

using std::size_t;

/******************************************************************************/

struct indval{
  int i;
  double x;
};

/******************************************************************************/

// Read/write memory-mapping
class SFBM {
public:
  SFBM(std::string path, size_t n, size_t m, std::vector<int> p) : n(n), m(m), p(p) {

    std::error_code error;
    this->ro_mmap.map(path, error);
    if (error) Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());

    this->data = reinterpret_cast<const indval*>(ro_mmap.data());
  }

  template<class C>
  C prod(const C& x) {

    C res(n);
    for (int i = 0; i < n; i++) res[i] = 0;

    for (int j = 0; j < m; j++) {

      auto lo = data + p[j];
      auto up = data + p[j + 1];

      for (auto it = lo; it < up; it++) {
        res[it->i] += it->x * x[j];
      }
    }

    return res;
  }

  template<class C>
  C crossprod(const C& x) {

    C res(m);

    for (int j = 0; j < m; j++)
      res[j] = dot_col(j, x);

    return res;
  }

  template<class C>
  double dot_col(int j, const C& x) {

    auto lo = data + p[j];
    auto up = data + p[j + 1];

    double cp = 0;
    for (auto it = lo; it < up; it++) {
      cp += it->x * x[it->i];
    }

    return cp;
  }

private:
  mio::mmap_source ro_mmap;
  const indval * data;
  size_t n;
  size_t m;
  std::vector<int> p;
};

/******************************************************************************/

#endif // SFBM_H
