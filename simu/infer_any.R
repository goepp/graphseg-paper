suppressMessages({library(tidyverse); library(gdata); library(RColorBrewer)})
suppressMessages({library(rgeos); library(rgdal)})
suppressMessages({library(glmnet); library(igraph); library(sf); library(Matrix)})
suppressMessages({library(pryr); library(limSolve)}) #library(graphseg)
suppressMessages({library(scales); library(latex2exp); library(flsa)})
suppressMessages({library(kableExtra); library(parallel); library(reshape2); library(dendextend)})
# library(RANN)
source("sf2nb.R")
source("div_pal.R")
source("infer_functions.R")
width_pdf = 5
height_pdf = 5
pal_div <- (c(brewer.pal(n = 9, name = "Blues") %>% rev,
              brewer.pal(n = 9, name = "Oranges"))) %>% colorRampPalette

# load data ---------------------------------------------------------------
args <- commandArgs(TRUE)
region <- ifelse(!is.na(args[1]), args[1], "utrecht")
area <- ifelse(!is.na(args[2]), args[2], "neigh")
zone <- ifelse(!is.na(args[3]), args[3], "district")
zone_code <- paste0(zone, "_code")

dataset <- readRDS(file = paste0("synthetic/", region, "_", area, ".rds"))
graph <- readRDS(file = paste0("synthetic/graph_", region, "_", area, ".rds"))

set.seed(87)
baseline <- data.frame("zone" = unique(dataset[[zone_code]]),
                       "baseline" = rpois(length(unique(dataset[[zone_code]])), 10))
dataset <- left_join(baseline, dataset, by = c("zone" = zone_code))

# Generate noisy signal ---------------------------------------------------------------
nrep <- 100
sigma_sq <- c(seq(0.1, 1.2, 0.1), 2, 3, 5)
names(sigma_sq) <- sigma_sq
set.seed(92)
gamma_array <- replicate(
  nrep, sapply(sigma_sq, function(sigma_sq) dataset$baseline + rnorm(nrow(dataset), 0, sigma_sq))
)

## Infer_res: computer rms and model_dim for each criterion ------
lambda <- 10 ^ seq(-3, 3, length = 50)

### Define functions
agraph <- function(gamma, graph, lambda = 10 ^ seq(-4, 4, length.out = 50),
                   weights = NULL, shrinkage = TRUE,
                   delta = 1e-10, tol = 1e-8,
                   thresh = 0.01, itermax = 50000) {
  p <- length(gamma)
  if (is.null(weights)) {
    weights <- rep(1, p)
  }
  precision <- Matrix::Diagonal(x = weights)
  nll <- model_dim <- rep(0, length(lambda))
  result <- matrix(NA, length(lambda), p)
  ind <- 1
  iter <- 1
  if (!(class(graph) %in% "igraph")) {
    stop("Graph must be an igraph object.")
  }
  if (length(gamma) != igraph::vcount(graph)) {
    stop("gamma must be a vector of length the number of vertices in the graph.")
  }
  edgelist_tmp <- igraph::as_edgelist(graph, names = FALSE)
  edgelist <- edgelist_tmp[order(edgelist_tmp[, 2]), c(2, 1)]
  adj <- as(igraph::as_adjacency_matrix(graph), "symmetricMatrix")
  sel <- adj
  converge <- FALSE
  weighted_laplacian_init <- lambda[ind] * (Matrix::Diagonal(x = Matrix::colSums(adj)) - adj) +
    Matrix::Diagonal(x = weights)
  chol_init <- Matrix::Cholesky(weighted_laplacian_init)
  while (iter < itermax) {
    sel_old <- sel
    weighted_laplacian <- lambda[ind] * (Diagonal(x = Matrix::colSums(adj)) - adj) +
      Matrix::Diagonal(x = weights)
    chol <- Matrix::update(chol_init, weighted_laplacian)
    theta <- Matrix::solve(chol, weights * gamma)
    adj@x <- 1 / ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    sel@x <- (theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 /
      ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    converge <- all(abs(sel@x - sel_old@x) < tol)
    if (converge) {
      graph_del <- igraph::delete_edges(graph, which((sel@x > 1 - thresh)[order(edgelist[, 2])]))
      segmentation <- igraph::components(graph_del)$membership
      if (shrinkage) {
        result[ind, ] <- stats::ave(as.vector(theta), segmentation)
      } else {
        result[ind, ] <- stats::ave(as.vector(gamma), segmentation)
      }
      nll[ind] <- 1 / 2 * t(result[ind, ] - gamma) %*% precision %*% (result[ind, ] - gamma)
      model_dim[ind] <- sum(diag(Matrix::solve(weighted_laplacian, precision)))
      ind <- ind + 1
    }
    iter <- iter + 1
    if (ind > length(lambda)) break
    # if (model_dim == 1) break
  }
  bic <- 2 * nll + log(p) * model_dim
  aic <- 2 * nll + 2 * model_dim
  gcv <- 2 * nll / (p * (1 - model_dim / p) ^ 2)
  return(list(result = result, bic = bic, gcv = gcv, model_dim = model_dim,
              nll = nll, aic = aic))
}
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
cl <- makeCluster(getOption("cl.cores", detectCores() - 2))
parallel::clusterExport(cl = cl, varlist = ls())
parallel::clusterEvalQ(cl = cl, {
  library(tidyverse); library(gdata)
  library(rgeos); library(rgdal)
  library(glmnet); library(igraph); library(sf); library(Matrix)
  library(pryr); library(limSolve); library(flsa)})
## Perform simulation ------
# Aridge
infer_res <- parApply(cl = cl, gamma_array, c(2, 3), infer, 
                      graph = graph, lambda = lambda,
                      baseline = dataset$baseline)
# Flsa
infer_res_flsa <- parApply(cl = cl, gamma_array, c(2, 3), infer_flsa, 
                      graph = graph, lambda = lambda,
                      baseline = dataset$baseline)
stopCluster(cl = cl)

### Mise en forme 
rms <- apply(infer_res, c(1, 2), function(z) unlist(z[[1]]$rms)) %>%
  "colnames<-"(sigma_sq)
rms_flsa <- apply(infer_res_flsa, c(1, 2), function(z) unlist(z[[1]]$rms)) %>%
  "colnames<-"(sigma_sq)
rms_df <- bind_rows(
  reshape2::melt(rms) %>%
    transmute("criterion" = Var1, "sigma_sq" = as.factor(Var2), "ind_rep" = Var3, "rms" = value) %>%
    mutate("treatment" = "aridge"),
  reshape2::melt(rms_flsa) %>%
    transmute("criterion" = Var1, "sigma_sq" = as.factor(Var2), "ind_rep" = Var3, "rms" = value) %>%
    mutate("treatment" = "flsa")
)
### Dimension
dim <- apply(infer_res, c(1, 2), function(z) unlist(z[[1]]$dim)) %>%
  "colnames<-"(sigma_sq)
dim_flsa <- apply(infer_res_flsa, c(1, 2), function(z) unlist(z[[1]]$dim)) %>%
  "colnames<-"(sigma_sq)
dim_df <- bind_rows(
  reshape2::melt(dim) %>%
    transmute("criterion" = Var1, "sigma_sq" = as.factor(Var2), "ind_rep" = Var3, "dim" = value) %>%
    mutate("treatment" = "aridge", "criterion" = criterion),
  reshape2::melt(dim_flsa) %>%
    transmute("criterion" = Var1, "sigma_sq" = as.factor(Var2), "ind_rep" = Var3, "dim" = value) %>%
    mutate("treatment" = "flsa"))

dim_rms_df <- inner_join(dim_df, rms_df, by = c("sigma_sq", "criterion", "ind_rep", "treatment"))
saveRDS(dim_rms_df, file = paste0("simu_results/cluster/", region, "_", area, "_pc_", zone, "_dim_rms_df.rds"))
