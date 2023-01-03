spmat2 <- as(cor(iris[1:4]), "dsCMatrix")
(X2 <- as_SFBM_corr_compact(spmat2))
(bin <- readBin(X2$sbk, what = integer(), size = 2, n = 100))
matrix(bin / 32767, 4)
spmat2

sp_prodVec(X2, rep(1, 4))
sp_cprodVec(X2, rep(1, 4))
Matrix::colSums(spmat2)

sp_solve_sym(X2, rep(1, 4))
X2[]
X2[, 1]
X2[1, ]
