svd_result <- svd(X.src[[1]])

U <- svd_result$u
d <- svd_result$d
V <- svd_result$v
m <- nrow(X.src[[1]])
n <- ncol(X.src[[1]])

# Compression
U[,80:min(n,m)] <- 0
d[80:min(n,m)] <- 0

Sigma <- diag(d)

V[,80:min(n,m)] <- 0

X_new <- Sigma%*%t(V)
y_new <- ginv(Sigma)%*%t(V) %*% t(X.src[[1]]) %*% y.src[[1]]
