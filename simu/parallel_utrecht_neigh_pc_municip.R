## NOTE: this script is present for legacy and is not used for the paper as of 
## May 2021

library(tidyverse); library(gdata); library(RColorBrewer)
library(rgeos); library(rgdal)
library(glmnet); library(igraph); library(sf); library(Matrix)
library(pryr); library(limSolve); library(graphseg)
library(scales); library(latex2exp); library(flsa)
library(kableExtra); library(parallel)

source("sf2nb.R")
source("div_pal.R")
source("infer_functions.R")
width_pdf = 5
height_pdf = 5
pal_div <- (c(brewer.pal(n = 9, name = "Blues") %>% rev,
              brewer.pal(n = 9, name = "Oranges"))) %>% colorRampPalette

load("synthetic/utrecht_neigh.RData")

set.seed(87)
baseline <- data.frame("municip_code" = levels(utrecht_neigh$municip_code),
                       "baseline" = rpois(nlevels(utrecht_neigh$municip_code), 10))
utrecht_neigh <- utrecht_neigh %>% left_join(baseline, by = "municip_code")

agraph <- function(gamma, graph, lambda = 10 ^ seq(-10, 2, length.out = 100),
                   weights = NULL, shrinkage = FALSE,
                   delta = 1e-10, tol = 1e-8,
                   thresh = 0.01, itermax = 50000) {
  gamma <- as.vector(gamma)
  lambda <- as.vector(lambda)
  p <- length(gamma)
  if (is.null(weights)) {
    weights <- rep(1, p)
  }
  bic <- rep(0, length(lambda))
  prec <- Diagonal(x = weights)
  nll <- model_dim <- gcv <- aic <- rep(0, length(lambda))
  result <- matrix(NA, length(lambda), p)
  ind <- 1
  iter <- 1
  if (!(class(graph) %in% "igraph")) {
    stop("Graph must be in igraph format")
  }
  edgelist_tmp <- as_edgelist(graph, names = FALSE)
  edgelist <- edgelist_tmp[order(edgelist_tmp[, 2]), c(2, 1)]
  adj <- as(as_adjacency_matrix(graph), "symmetricMatrix")
  sel <- adj
  converge <- FALSE
  weighted_laplacian_init <- lambda[ind] * (Diagonal(x = colSums(adj)) - adj) + Diagonal(x = weights)
  chol_init <- Cholesky(weighted_laplacian_init)
  while (iter < itermax) {
    sel_old <- sel
    weighted_laplacian <- lambda[ind] * (Diagonal(x = colSums(adj)) - adj) + Diagonal(x = weights)
    chol <- update(chol_init, weighted_laplacian)
    theta <- solve(chol, weights * gamma)
    adj@x <- 1 / ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    sel@x <- (theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 /
      ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    converge <- all(abs(sel@x - sel_old@x) < tol)
    if (converge) {
      graph_del <- delete_edges(graph, which((sel@x > 1 - thresh)[order(edgelist[, 2])]))
      segmentation <- components(graph_del)$membership
      if (shrinkage) {
        result[ind, ] <- ave(as.vector(theta), segmentation)
      } else {
        result[ind, ] <- ave(as.vector(gamma), segmentation)
      }
      nll[ind] <- 1 / 2 * t(result[ind, ] - gamma) %*% prec %*% (result[ind, ] - gamma)
      model_dim[ind] <- sum(diag(solve(weighted_laplacian, prec)))
      bic[ind] <- 2 * nll[ind] + log(p) * model_dim[ind]
      aic[ind] <- 2 * nll[ind] + 2 * model_dim[ind]
      gcv[ind] <- 2 * nll[ind] / (p * (1 - model_dim[ind] / p) ^ 2)
      ind <- ind + 1
    }
    iter <- iter + 1
    if (ind > length(lambda)) break
    # if (model_dim == 1) break
  }
  return(list(result = result, bic = bic, gcv = gcv, model_dim = model_dim,
              nll = nll, aic = aic))
}
graph2connlist <- function(graph) {
  connlist <- igraph::as_adj_list(graph) %>%
    lapply(as.character) %>%
    lapply(as.numeric) %>%
    lapply(function(a) a - 1) %>%
    lapply(function(a) {
      attributes(a) <- NULL
      return(a)
    }) %>%
    lapply(as.integer) %>%
    'class<-'("connListObj")
  return(connlist)
}
#' @export
flsa_graph <- function(gamma, graph, lambda) {
  flsa_fit <- flsa::flsa(gamma, connListObj = graph2connlist(graph),
                         lambda2 = lambda)
  model_dim <- apply(flsa_fit, 1, function(a) length(unique(a))) %>% unname
  nll <- 1 / 2 * apply(flsa_fit - t(replicate(length(lambda), gamma)),
                       1,
                       function(a) sum(a ^ 2))
  list("result" = flsa_fit, "aic" = 2 * nll + 2 * model_dim,
       "bic" = 2 * nll + log(length(gamma)) * model_dim,
       "gcv" = 2 * nll / ((length(gamma) - model_dim) ^ 2),
       "model_dim" = model_dim, "nll" = nll)
}

## Simulate data
nrep <- 100
# sigma_sq <- exp(seq(log(0.05), log(20), length.out = 10))
sigma_sq <- c(seq(0.1, 1.2, 0.1), 2, 3, 5)
names(sigma_sq) <- sigma_sq

set.seed(92)
gamma_array <- replicate(
  nrep, sapply(sigma_sq, function(sigma_sq) utrecht_neigh$baseline + rnorm(nrow(utrecht_neigh), 0, sigma_sq))
)

## Perform simulation
lambda <- 10 ^ seq(-4, 3, length = 50)

cl <- makeCluster(getOption("cl.cores", 10))
parallel::clusterExport(cl = cl, varlist = ls())
parallel::clusterEvalQ(cl = cl, {
             library(tidyverse); library(gdata); library(RColorBrewer)
             library(rgeos); library(rgdal)
             library(glmnet); library(igraph); library(sf); library(Matrix)
             library(pryr); library(limSolve)
             library(scales); library(latex2exp); library(flsa)
             library(kableExtra)})
infer_res_flsa2 <- parApply(cl, gamma_array[, , 10:20], c(2, 3), infer_flsa, graph = graph_utrecht_neigh,
                           lambda = lambda, baseline = utrecht_neigh$baseline)
#1:10 infer_res_flsa
#10:20 infer_res_flsa2
stopCluster(cl = cl)

dim(infer_res_flsa)
# infer_res_flsa <- parApply(cl, gamma_array[, , 1:2], c(2, 3), infer_flsa, graph = graph_utrecht_neigh,
                           # lambda = lambda, baseline = utrecht_neigh$baseline)

parApply(cl, gamma_array[, 1:2, 1:2], c(2, 3), flsa_graph, graph = graph_utrecht_neigh,
         lambda = lambda)
