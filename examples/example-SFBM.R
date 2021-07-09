spmat2 <- Matrix::Diagonal(4, 0:3)
spmat2[4, 2] <- 5
spmat2[1, 4] <- 6
spmat2[3, 4] <- 7
spmat2

# Stores all (i, x) for x != 0
(X2 <- as_SFBM(spmat2))
matrix(readBin(X2$sbk, what = double(), n = 100), 2)

# Stores only x, but all (even the zero ones) from first to last being not 0
(X3 <- as_SFBM(spmat2, compact = TRUE))
X3$first_i
readBin(X3$sbk, what = double(), n = 100)
