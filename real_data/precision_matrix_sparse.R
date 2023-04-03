library(ggplot2)
library(magrittr)
library(glasso)
library(Matrix)
library(huge)

sigma <- scan(file = "V_mrf.txt") %>% matrix(ncol = sqrt(length(.))) / 1e8
isSymmetric(sigma)
image(sigma[1:100, 1:100])

## Check the high covariance values
mean(diag(sigma)) # Very low variance
corr <- sigma[col(sigma) > row(sigma)]
corr <- corr[corr != 0]
median(abs(corr)) # Quite low on average


## Check that the diagonal elements of sigma concur
## with the standard errors previously used
se <- sqrt(diag(sigma))
mrf <- read.delim(file = "mrfterm_overweight_utrecht.txt")
sigma_sq <- mrf$se.mrf
all.equal(sigma_sq, se, tol = 1e-6)


## Correlation
cor <- cov2cor(sigma)
cor_dp <- cor + 20 * Diagonal(n = ncol(cor))
# small <- cor_dp[1:20, 1:20] %>% as.matrix()

sigma_dp <- sigma + 0.01 * Diagonal(n = ncol(sigma))
all(eigen(sigma_dp[1:50, 1:50])$values > 0)

## Compare two packages for graphical lasso
## on a subset of sigma
small <- sigma[1:100, 1:100] %>% as.matrix()
lambda <- 0.005

res_glasso_huge <- huge.glasso(small, lambda = lambda)
prec_1 <- res_glasso_huge$icov[[1]] %>% Matrix(sparse = TRUE)
image(prec_1)

res_glasso <- glasso(small, rho = lambda)
prec_2 <- res_glasso$wi %>% Matrix(sparse = TRUE)
image(prec_2)

all.equal(res_glasso_huge$icov[[1]],
          res_glasso$wi)
image(Matrix(res_glasso_huge$icov[[1]] - res_glasso$wi, sparse = TRUE))

## Verify that huge.glasso computes the inverse of
## the covariance matrix
small <- sigma[1:100, 1:100] %>% as.matrix()
lambda <- 0.01

res_glasso_huge <- huge.glasso(small, lambda = lambda)
prec_1 <- res_glasso_huge$icov[[1]] %>% Matrix(sparse = TRUE)
image(prec_1)

res_glasso <- glasso(small, rho = lambda)
prec_2 <- res_glasso$wi %>% Matrix(sparse = TRUE)
image(prec_2)

plot(diag(prec_1), 1 / diag(small))

## Compute the precision matrix
small <- sigma[1:2000, 1:2000] %>% as.matrix()
lambda <- 0.0001

tmp <- huge.glasso(sigma, lambda = lambda, scr = TRUE)
prec <- tmp$icov[[1]] %>% Matrix(sparse = TRUE)
image(prec)
prec <- prec %>% forceSymmetric()
# save(prec, file = "utrecht_prec.RData")

all.equal(as.matrix(prec %*% sigma), diag(ncol(sigma)))
norm(as.matrix(prec %*% sigma) - diag(ncol(sigma))) / norm(as.matrix(prec %*% sigma))
