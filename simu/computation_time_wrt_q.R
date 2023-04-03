suppressMessages({library(tidyverse); library(gdata); library(RColorBrewer)})
suppressMessages({library(rgeos); library(rgdal); library(microbenchmark)})
suppressMessages({library(glmnet); library(igraph); library(sf); library(Matrix)})
suppressMessages({library(pryr); library(limSolve)})
##library(graphseg)
suppressMessages({library(scales); library(latex2exp); library(flsa)})
suppressMessages({library(kableExtra); library(microbenchmark); library(parallel)})
#library(genlasso)

## Load utils functions
source(file = "utils/infer_functions.R")

## Method with automatically extracted subgraphs -----
load(file = "simu/synthetic/graph_sub_q.RData")
lambda <- 10 ^ seq(-4, 4, length = 50)
n_times <- 100
bench_agraph_flsa <- function(signal, graph, lambda, times) {
  microbenchmark(
    "agraph" = agraph(signal, graph, lambda),
    "flsa" = flsa_graph(signal, graph, lambda), times = times) %>% as.data.frame()
}

print(paste("NUMBER OF CORES:", detectCores()))


# Do benchmark ------------------------------------------------------------
## Option 1: do benchmark on local computer
# bench_df_q <- lapply(list(sub_q_10$signal, sub_q_30$signal,
#                         sub_q_100$signal, sub_q_300$signal),
#                    bench_agraph_flsa,
#                    graph = graph_sub_q,
#                    lambda = lambda, time = n_times) %>% do.call(rbind, .)

# saveRDS(bench_df_q, file = "simu/simu_results/bench_df_q.rds")

# ## Option 2: do benchmark on cluster
# cl <- makeCluster(getOption("cl.cores", detectCores()))
# parallel::clusterExport(cl = cl, varlist = ls())
# parallel::clusterEvalQ(cl = cl, {
#   library(tidyverse)
#   library(rgeos); library(rgdal)
#   library(glmnet); library(igraph); library(sf); library(Matrix)
#   library(pryr); library(limSolve); library(microbenchmark)
#   library(flsa); 
#   agraph <- function(gamma, graph, lambda = 10 ^ seq(-4, 4, length.out = 50),
#                      weights = NULL, shrinkage = TRUE,
#                      delta = 1e-10, tol = 1e-8,
#                      thresh = 0.01, itermax = 50000) {
#     p <- length(gamma)
#     if (is.null(weights)) {
#       weights <- rep(1, p)
#     }
#     precision <- Matrix::Diagonal(x = weights)
#     nll <- model_dim <- rep(0, length(lambda))
#     result <- matrix(NA, length(lambda), p)
#     ind <- 1
#     iter <- 1
#     if (!(class(graph) %in% "igraph")) {
#       stop("Graph must be an igraph object.")
#     }
#     if (length(gamma) != igraph::vcount(graph)) {
#       stop("gamma must be a vector of length the number of vertices in the graph.")
#     }
#     edgelist_tmp <- igraph::as_edgelist(graph, names = FALSE)
#     edgelist <- edgelist_tmp[order(edgelist_tmp[, 2]), c(2, 1)]
#     adj <- as(igraph::as_adjacency_matrix(graph), "symmetricMatrix")
#     sel <- adj
#     converge <- FALSE
#     weighted_laplacian_init <- lambda[ind] * (Matrix::Diagonal(x = Matrix::colSums(adj)) - adj) +
#       Matrix::Diagonal(x = weights)
#     chol_init <- Matrix::Cholesky(weighted_laplacian_init)
#     while (iter < itermax) {
#       sel_old <- sel
#       weighted_laplacian <- lambda[ind] * (Diagonal(x = Matrix::colSums(adj)) - adj) +
#         Matrix::Diagonal(x = weights)
#       chol <- Matrix::update(chol_init, weighted_laplacian)
#       theta <- Matrix::solve(chol, weights * gamma)
#       adj@x <- 1 / ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
#       sel@x <- (theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 /
#         ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
#       converge <- all(abs(sel@x - sel_old@x) < tol)
#       if (converge) {
#         graph_del <- igraph::delete_edges(graph, which((sel@x > 1 - thresh)[order(edgelist[, 2])]))
#         segmentation <- igraph::components(graph_del)$membership
#         if (shrinkage) {
#           result[ind, ] <- stats::ave(as.vector(theta), segmentation)
#         } else {
#           result[ind, ] <- stats::ave(as.vector(gamma), segmentation)
#         }
#         nll[ind] <- 1 / 2 * t(result[ind, ] - gamma) %*% precision %*% (result[ind, ] - gamma)
#         model_dim[ind] <- sum(diag(Matrix::solve(weighted_laplacian, precision)))
#         ind <- ind + 1
#       }
#       iter <- iter + 1
#       if (ind > length(lambda)) break
#       # if (model_dim == 1) break
#     }
#     bic <- 2 * nll + log(p) * model_dim
#     aic <- 2 * nll + 2 * model_dim
#     gcv <- 2 * nll / (p * (1 - model_dim / p) ^ 2)
#     return(list(result = result, bic = bic, gcv = gcv, model_dim = model_dim,
#                 nll = nll, aic = aic))
#   };
#   flsa_graph <- function(gamma, graph, lambda) {
#     flsa_fit <- flsa::flsa(gamma, connListObj = graph2connlist(graph),
#                            lambda2 = lambda)
#     model_dim <- apply(flsa_fit, 1, function(a) length(unique(a))) %>% unname
#     nll <- 1 / 2 * apply(flsa_fit - t(replicate(length(lambda), gamma)),
#                          1,
#                          function(a) sum(a ^ 2))
#     list("result" = flsa_fit, "aic" = 2 * nll + 2 * model_dim,
#          "bic" = 2 * nll + log(length(gamma)) * model_dim,
#          "gcv" = 2 * nll / ((length(gamma) - model_dim) ^ 2),
#          "model_dim" = model_dim, "nll" = nll)
#   };
#   graph2connlist <- function(graph) {
#     connlist <- igraph::as_adj_list(graph) %>%
#       lapply(as.character) %>%
#       lapply(as.numeric) %>%
#       lapply(function(a) a - 1) %>%
#       lapply(function(a) {
#         attributes(a) <- NULL
#         return(a)
#       }) %>%
#       lapply(as.integer) %>%
#       'class<-'("connListObj")
#     return(connlist)
#   }})
# set.seed(45)
# bench_df <- clusterMap(cl = cl, bench_agraph_flsa,
#                        list(sub_10000$signal),
#                        list(graph_sub_10000),
#                        MoreArgs = list(lambda = lambda, times = n_times),
#                        SIMPLIFY = FALSE) %>%
#   do.call(rbind, .)
# stopCluster(cl)
# saveRDS(bench_df, file = "simu_results/bench_df3.rds")

# bench_df <- readRDS(file = "simu/simu_results/bench_df.rds") # OLD
# bench_df <- readRDS(file = "simu/simu_results/bench_df_cluster.rds")

# create runtime table -----------------------------------------------------
# bench_df_q <- readRDS(file = "simu/simu_results/bench_df_cluster.rds")
temp <- bench_df_q %>% 
  mutate(sample_size = c(rep(c(10, 30, 100, 300), each = 2 * 100)))
bench_margin_df <- temp %>%
  dplyr::group_by(expr, sample_size) %>%
  dplyr::summarize(mean = mean(time),
                   sd = sd(time)) %>%
  mutate(mean = sprintf("%.2f", signif(mean, digits = 2)),
         sd = sprintf("%.2f", signif(sd, digits = 2))) %>%
  mutate(value_err = paste('$', mean, "\\pm", sd, '$'))
min_col <- bench_margin_df %>% select(-c(value_err, sd)) %>%
  reshape2::acast(sample_size ~ expr) %>%
  apply(1, function(a) which.min(a))
prod_table <- bench_margin_df %>% reshape2::acast(sample_size ~ expr) %>%
  as.data.frame(stringsAsFactors = FALSE)
for (ind in seq_along(min_col)) {
  prod_table[ind, min_col[ind]] <- paste0("\\boldmath{",
                                          prod_table[ind, min_col[ind]],
                                          "}")
}
# bench_kable <- prod_table %>%
#   tibble::rownames_to_column(var = "sample_size") %>%
#   kable("latex", booktabs = TRUE, escape = FALSE,
#         caption = "Computing time in seconds of both method for datasets of different sample sizes.",
#         col.names = linebreak(c("Sample \n size",
#                                 "Adaptive \n ridge", "FLSA"), align = "c")
#   ) %>%
#   kable_styling(latex_options = c("striped"), position = "center")
# 
# bench_kable %>% cat(file = "table/computation_time.tex")

# create runtime graph ----------------------------------------------------
bench_df_q <- readRDS(file = "simu/simu_results/bench_df_q.rds")
# temp <- bench_df_q %>%
#   mutate(sample_size = rep(c(10, 30, 100, 300), each = 2 * n_times))
temp <- bench_df_q %>% 
  mutate(q = c(rep(c(10, 30, 100, 300), each = 2 * 100)))
bench_df_for_graph <- temp %>% 
  mutate(time = time / 1e9) %>% 
  dplyr::group_by(expr, q)

ggplot(bench_df_for_graph, aes(as.factor(q), time, color = expr)) +
  geom_boxplot() +
  theme_minimal() +
  scale_color_discrete(labels = c("adaptive ridge", "flsa")) +
  ylab("Computing time (seconds)") + 
  xlab(TeX("Number of zones (q)")) + 
  ggtitle("With p = 3000 areas") +
  # scale_y_log10() +
  theme(legend.position = c(0.2, 0.2), 
        text = element_text(family = "LM Roman 10", face = "plain", size = 14),
        legend.title = element_blank())
width_pdf <- 5
height_pdf <- 5

# ggsave(file = "simu/figure/computing_time_boxplot_linear_scale.pdf", 
#        device = cairo_pdf,
#        width = width_pdf, height = height_pdf)

### With crossbar instead of boxplots (allows for numeric x axis)
bench_df_crossbar <- bench_df_for_graph %>% 
  group_by(q, expr) %>% 
  summarize(mean = mean(time), std = sd(time), .groups = "drop") %>% 
  mutate(lower = mean - std, upper = mean + std)

ggplot(bench_df_crossbar, aes(x = q, y = mean, ymin = lower,
                              ymax = upper, color = expr)) +
  geom_crossbar() + geom_line() +
  theme_minimal() +
  scale_color_discrete(labels = c("adaptive ridge", "flsa")) +
  ylab("Computing time (seconds)") + 
  xlab(TeX("Number of zones (q)")) +  
  ggtitle("With p = 3000 areas") +
  # scale_y_log10() +
  theme(legend.position = c(0.8, 0.9), 
        text = element_text(family = "LM Roman 10", face = "plain", size = 14),
        legend.title = element_blank())
width_pdf <- 6
height_pdf <- 4

ggsave(file = "simu/figure/computing_time_q.pdf",
       device = cairo_pdf,
       width = width_pdf, height = height_pdf)


