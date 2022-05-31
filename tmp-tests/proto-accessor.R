corr0 <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))
library(bigsparser)

corr <- as_SFBM(corr0)
sbk <- corr$backingfile
bk <- sub("\\.sbk$", ".bk", sbk)
file.copy(sbk, bk)
library(bigstatsr)
X <- FBM(nrow = 2, ncol = corr$nval, backingfile = sub_bk(bk),
         create_bk = FALSE, is_read_only = TRUE)
X[, 1:7]
corr0[1:10, 1:2, drop = FALSE]

Rcpp::sourceCpp('src/spmat_accessor.cpp')
ind_row <- 1:10
ind_col <- 1:2
res <- access_subset(corr, ind_row, ind_col)

sub <- new("dgCMatrix")
sub@Dim <- c(length(ind_row), length(ind_col))
sub@i <- res$i
sub@p <- as.integer(res$p)
sub@x <- res$x
sub
