library(bigstatsr)
library(testthat)
spmat0 <- Matrix::rsparsematrix(1000, 1000, 0.01, symmetric = TRUE)
X <- as_SFBM(spmat0)
spmat <- as(spmat0, "dgCMatrix")
X2 <- as_SFBM(spmat)
expect_identical(readBin(X$sbk,  what = 1, n = 1e6),
                 readBin(X2$sbk, what = 1, n = 1e6))
expect_equal(X$p, spmat@p)

col_count <- bigsparser:::col_count_sym(spmat0@p, spmat0@i)
head(col_count)
sum(col_count)
