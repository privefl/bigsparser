## bigsparser: sparse matrix format with data on disk

### Features

For now, only a few features are implemented:

- convert a dgCMatrix to an SFBM, a Sparse Filebacked Big Matrix

- compute the product and crossproduct of an SFBM with a vector

- solve Ax=b, where A is a symmetric SFBM and b is a vector

### Installation

```r
remotes::install_github("privefl/bigsparser")
```
