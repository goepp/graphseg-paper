infer_any <- function(gamma, method, graph, lambda, baseline, true_labels, 
                      score = mclust::adjustedRandIndex) {
  if (!(method %in% c("agraph", "flsa"))) {
    stop("Error: method must be agraph of flsa")
  }
  if (method == "agraph") {
    theta <- agraph(gamma = gamma,
                    graph = graph,
                    lambda = lambda,
                    weights = NULL,
                    shrinkage = TRUE,
                    tol = 1e-8,
                    itermax = 50000)
  } else if (method == "flsa") {
    theta <- flsa_graph(gamma = gamma,
                        graph = graph,
                        lambda = lambda)
  }
  # RMS
  rms_aic <- sqrt(sum((theta$result[which.min(theta$aic), ] - baseline) ^ 2) / length(gamma))
  rms_bic <- sqrt(sum((theta$result[which.min(theta$bic), ] - baseline) ^ 2) / length(gamma))
  rms_gcv <- sqrt(sum((theta$result[which.min(theta$gcv), ] - baseline) ^ 2) / length(gamma))
  # Model dimension
  dim_aic <- theta$model_dim[which.min(theta$aic)]
  dim_bic <- theta$model_dim[which.min(theta$bic)]
  dim_gcv <- theta$model_dim[which.min(theta$gcv)]
  # Clustering score
  clust_score_aic <- score(as.numeric(as.factor(theta$result[which.min(theta$aic), ])), true_labels)
  clust_score_bic <- score(as.numeric(as.factor(theta$result[which.min(theta$bic), ])), true_labels)
  clust_score_gcv <- score(as.numeric(as.factor(theta$result[which.min(theta$gcv), ])), true_labels)
  return(list("rms" = c("aic" = rms_aic, "bic" = rms_bic, "gcv" = rms_gcv),
              "dim" = c("aic" = dim_aic, "bic" = dim_bic, "gcv" = dim_gcv),
              "clust_score" = c("aic" = clust_score_aic, "bic" = clust_score_bic, "gcv" = clust_score_gcv)))
}
infer <- function(gamma, graph, lambda, baseline, true_labels, 
                  score = mclust::adjustedRandIndex) {
  infer_any(gamma, method = "agraph", graph, lambda, baseline, true_labels, 
            score = mclust::adjustedRandIndex)
}
infer_flsa <- function(gamma, graph, lambda, baseline, true_labels, 
                       score = mclust::adjustedRandIndex) {
  infer_any(gamma, method = "flsa", graph, lambda, baseline, true_labels, 
            score = mclust::adjustedRandIndex)
}
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
agraph_prec <- function(gamma, graph, prec,
                        lambda = 10 ^ seq(-4, 4, length.out = 50),
                        weights = NULL, shrinkage = TRUE,
                        delta = 1e-10, tol = 1e-8,
                        thresh = 0.01, itermax = 10000,
                        average_fast = FALSE) {
  gamma <- as.vector(gamma)
  lambda <- as.vector(lambda)
  if (!isSymmetric(prec)) {
    prec <- Matrix::forceSymmetric(prec)
    warning("Precision matrix was not symmetric. It was coerced to symmetric.")
  }
  if (!"sparseMatrix" %in% is(prec)) { # prec is not sparse
    prec <- as(prec, "symmetricMatrix")
  }
  p <- length(gamma)
  prec_gamma <- prec %*% gamma
  nll <- model_dim <- rep(0, length(lambda))
  result <- matrix(NA, length(lambda), p)
  ind <- 1
  iter <- 1
  if (!(class(graph) %in% "igraph")) {
    stop("Graph must be in igraph format")
  }
  edgelist_tmp <- igraph::as_edgelist(graph, names = FALSE)
  edgelist <- edgelist_tmp[order(edgelist_tmp[, 2]), c(2, 1)]
  adj <- Matrix::forceSymmetric(igraph::as_adjacency_matrix(graph, sparse = TRUE))
  sel <- adj
  converge <- FALSE
  weighted_laplacian_init <- lambda[ind] * (Matrix::Diagonal(x = colSums(adj)) - adj) + prec
  chol_init <- Matrix::Cholesky(weighted_laplacian_init)
  while (iter < itermax) {
    sel_old <- sel
    weighted_laplacian <- lambda[ind] * (Matrix::Diagonal(x = colSums(adj)) - adj) + prec
    chol <- Matrix::update(chol_init, weighted_laplacian)
    theta <- Matrix::solve(chol, prec_gamma)
    adj@x <- 1 / ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    sel@x <- (theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 /
      ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    converge <- all(abs(sel@x - sel_old@x) < tol); converge
    if (converge) {
      graph_del <- igraph::delete_edges(graph, which((sel@x > 1 - thresh)[order(edgelist[, 2])]))
      segmentation <- igraph::components(graph_del)$membership
      if (shrinkage) {
        if (average_fast) {
          result[ind, ] <- stats::ave(as.vector(theta), segmentation)
        } else {
          result[ind, ] <- ave_matrix_weights(as.vector(theta), segmentation)
        }
      } else {
        result[ind, ] <- stats::ave(as.vector(gamma), segmentation)
      }
      nll[ind] <- 1 / 2 * t(result[ind, ] - gamma) %*% prec %*% (result[ind, ] - gamma)
      model_dim[ind] <- sum(diag(Matrix::solve(weighted_laplacian, prec)))
      ind <- ind + 1
    }
    iter <- iter + 1
    if (ind > length(lambda)) break
  }
  bic <- 2 * nll + log(p) * model_dim
  aic <- 2 * nll + 2 * model_dim
  gcv <- 2 * nll / ((p - model_dim) ^ 2)
  return(list(result = result, aic = aic, bic = bic,
              gcv = gcv, model_dim = model_dim, nll = nll))
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
weighted_mean <- function(vector_of_incides, x, prec) {
  weighted_sum = sum(prec[vector_of_incides, vector_of_incides] %*% x[vector_of_incides])
  factor = sum(prec[vector_of_incides, vector_of_incides])
  rep(weighted_sum / factor, length(vector_of_incides))
}
ave_matrix_weights <- function(x, levels) {
  split(seq_along(as.vector(theta)), levels) %>% 
    lapply(weighted_mean, x = as.vector(theta), prec = prec) %>% 
    unsplit(segmentation)
}
# Define clustering score as the Fowlkes-Mallows index
# score <- dendextend::FM_index
score <- mclust::adjustedRandIndex
