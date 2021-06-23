load(url("https://www.dropbox.com/s/c13uygnjh6yh7vf/to-test-ldpred2.RData?raw=1"))
dim(corr)
class(corr)  # dsCMatrix
class(corr[, 1:ncol(corr)])  # -> get a dgCMatrix -> could simply add blocks of columns


corr0 <- corr[, 1:ncol(corr)]

max(corr0@p)     # 2884089
length(corr0@i)  # same

corr0[, 1] <- 0
corr0 <- Matrix::drop0(corr0)

Rcpp::sourceCpp('tmp-tests/test-count-nonzero.cpp')
(counts <- col_count_compact(corr0@p, corr0@i))

sum(counts)
length(corr0@x)


corr2 <- Matrix::Diagonal(4, 0:3)
corr2[4, 2] <- 5
corr2[1, 4] <- 6
corr2[3, 4] <- 7
corr2

Rcpp::sourceCpp('tmp-tests/test-count-nonzero.cpp')
(counts2 <- col_count_compact(corr2@p, corr2@i))
p <- c(0, cumsum(counts2))
counts2
write_compact_val(tmp <- tempfile(), corr2@p, corr2@i, corr2@x)
readBin(tmp, what = 1, n = 100)

typeof(cumsum(seq_len(1e5)))  # /!\ overflow /!\ -> need +0

res <- col_range_sym(corr@p, corr@i)
ranges <- cbind(res[[1]], res[[2]])
plot(ranges); abline(0, 1, col = "red")

ranges_true <- t(apply(corr, 2, function(x) range(which(x != 0))))
plot(ranges_true)
identical(ranges + 1L, ranges_true)
