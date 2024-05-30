## bigsparser: sparse matrix format with data on disk

<!-- badges: start -->
[![R-CMD-check](https://github.com/privefl/bigsparser/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/privefl/bigsparser/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/bigsparser)](https://CRAN.R-project.org/package=bigsparser)
<!-- badges: end -->

### Features

For now, only a few features are implemented:

- convert a dgCMatrix or a dsCMatrix to an SFBM, a Sparse Filebacked Big Matrix

- grow an SFBM using `$add_columns()` (similar to 'cbind' or 'bdiag')

- compute the product and crossproduct of an SFBM with a vector

- solve Ax=b, where A is a symmetric SFBM and b is a vector

- access the subset of an SFBM as a dgCMatrix (matrix accessor, since v0.6.1)

A **new compact format** is available (since v0.5), which is useful when non-zero values in columns are contiguous (or almost).


### Installation

```r
# CRAN
install.packages("bigsparser")

# latest GitHub version
remotes::install_github("privefl/bigsparser")
```
