library(tidyverse); library(gdata); library(RColorBrewer)
library(rgeos); library(rgdal); library(lattice)
library(glmnet); library(igraph); library(sf)
library(mgcv); library(Matrix)
library(pryr); library(limSolve)
library(RANN); library(car)

source("utils/sf2nb.R")
source("utils/div_pal.R")
source("utils/infer_functions.R")

load(file = "real_data/utrecht_prec.RData")
width_pdf = 6
height_pdf = 6

# Load raw data -----------------------------------------------------------
mrf <- read.delim(file = "real_data/raw_data/mrfterm_overweight_utrecht.txt")
names(mrf)
gamma <- mrf$mrf
sigma_sq <- mrf$se.mrf
polygon <- st_read(dsn = "real_data/raw_data/map_netherland.geojson", quiet = TRUE) %>%
  filter(statcode %in% mrf$bu_code) %>%
  mutate(gamma = gamma)

# pal_seq <- brewer.pal(n = 9, name = "Blues") %>% colorRampPalette()
# pal_div <- diverging_hcl(50, palette = "Blue-Yellow 3")
# plot(polygon["gamma"], border = NA, nbreaks = 50, pal = pal_div, breaks = "fisher")

adj_sparse <- (st_intersects(polygon, polygon, sparse = FALSE) -
                 diag(nrow(polygon))) %>%
  Matrix(sparse = TRUE) %>%
  'dimnames<-'(list(polygon$statcode, polygon$statcode))
graph <- sf2nb(polygon, queen = FALSE, sparse = FALSE) %>%
  graph_from_adjacency_matrix(mode = "undirected")

# With covariance, with shrinkage ---------------------------------------------
# lambda <- 10 ^ seq(-4, -1, by = 0.03) # length = 101
lambda <- 10 ^ seq(-4, 2, by = 0.03) # length = 201

mrf_seg_prec <- agraph_prec(gamma = gamma,
                            graph = graph,
                            prec / 2, lambda = lambda, weights = NULL,
                            shrinkage = TRUE, tol = 1e-8, itermax = 50000,
                            average_fast = TRUE)
# saveRDS(mrf_seg_prec, file = "real_data/utrecht_mrf_seg_prec.rds")
# mrf_seg_prec <- readRDS(file = "real_data/utrecht_mrf_seg_prec.rds")

### Load the version of the paper sent to CSDA:
# mrf_seg_prec <- readRDS(file = "real_data/utrecht_mrf_seg_prec_save.rds")
 

# mrf_seg_prec$gcv %>% plot()
# mrf_seg_prec$gcv %>% log %>%  plot(x = log10(lambda))
# mrf_seg_prec$nll %>% plot()
# p <- nrow(polygon)
# plot(log(1 / (p * (1 - mrf_seg_prec$model_dim / p) ^ 2)))
# points(goodness_of_fit, col = "red")
# tt <- goodness_of_fit / (p * (1 - mrf_seg_prec$model_dim / p) ^ 2)
# plot(goodness_of_fit / (p * (1 - mrf_seg_prec$model_dim / p) ^ 2))
# 
# 
# mrf_seg_prec$result  %>% matplot(x = lambda, type = 'l', lwd = 0.5, log = "x")
# mrf_seg_prec$model_dim
# mrf_seg_prec$model_dim[which.min(mrf_seg_prec$gcv)]
# mrf_seg_prec$model_dim[which.min(mrf_seg_prec$aic)]
# mrf_seg_prec$model_dim[which.min(mrf_seg_prec$bic)]
# ind <- which.min(mrf_seg_prec$bic)
# ind <- 90
# ind <- which.min(mrf_seg_prec$aic)
# 
# polygon["seg_prec"] <- mrf_seg_prec$result[ind, ]
# col_seg_prec <- div_pal(polygon[["seg_prec"]], by = 0.01,
#                         lower = -0.3, upper = 0.4)$col
# pal_seg_prec <- div_pal(polygon[["seg_prec"]], by = 0.01,
#                         lower = -0.3, upper = 0.4)$pal
# (legend_breaks <- seq(-0.3, 0.4, length = 11)) # where the legend breaks come from
# legend_breaks_names <- mapply(function(a, b) paste0(a, " to ", b), legend_breaks[-1], legend_breaks[-11])
# plot(polygon["seg_prec"], col = col_seg_prec,
#      # main = paste("seg_prec:",
#      #              apply(mrf_seg_prec$result, 2, function(vec) length(unique(vec)))[ind],
#      #              "params"),
#      main = NA,
#      border = NA)
# legend("right",
#        legend = c("below -0.25", "-0.25 to -0.20", "-0.20 to -0.15", "-0.15 to -0.10",
#                   "-0.10 to -0.05", "-0.05 to 0.00", "0.00 to 0.05", "0.05 to 0.10",
#                   "0.10 to 0.15", "0.15 to 0.20", "0.20 to 0.25", "above 0.25"),
#        fill = pal_div(12),
#        bty = "n",
#        cex = 0.58,
#        ncol = 1)
# dev.copy2pdf(file = "utrecht_mrf_seg_prec.pdf")
# 
# ## Conclusion: aic, bic and gcv does not yield interesting estimates
# 
# # Plot unseg prevalence ---------------------------------------------------
# col_unseg = div_pal(polygon[["gamma"]], by = 0.01, lower = -0.3, upper = 0.4)$col
# breaks = div_pal(polygon[["gamma"]], by = 0.01, lower = -0.3, upper = 0.4)$breaks
# (legend_breaks <- seq(-0.3, 0.4, length = 11)) # where the legend breaks come from
# legend_breaks_names <- mapply(function(a, b) paste0(a, " to ", b), legend_breaks[-1], legend_breaks[-11])
# 
# # plot(polygon["gamma"], border = NA, col = col_unseg, main = NA)
# plot(polygon["gamma"], border = NA, col = col_unseg,
#      key.pos = get_key_pos(polygon["gamma"]),
#      key.length = 0.618 * 2,
#      key.width = lcm(1.8 * 2))
# legend("right",
#        legend = c("below -0.3", legend_breaks_names, "above 0.4"),
#        fill = pal_div(12), #bty = "n",
#        cex = 0.50, ncol = 1, border = NA)
# dev.copy2pdf(file = "real_data/figure/utrecht_mrf_unseg.pdf",
#              width = 5, height = 5)

# OLD -- Perform MRF estimation --------------------------------------------------
# lambda <- 10 ^ seq(-8, -4, length = 100)
# 
# mrf_res <- agraph(gamma,
#                   sigma_sq,
#                   adj_sparse, pen = lambda,
#                   maxiter = 10000)
# save(mrf_seg_slow, file = "mrf_seg_slow.RData")
# load(file = "mrf_seg_slow.RData")
# polygon["seg_slow"] <- mrf_seg_slow$par[[1]]
# plot(polygon["seg_slow"], nbreaks = 50, pal = pal_div, breaks = "fisher", lwd = 0.1)
# 
# mrf_seg_shrink <- agraph(gamma = gamma,
#                          graph = graph,
#                          lambda = lambda,
#                          weights = 1 / 2 * sigma_sq ^ 2,
#                          shrinkage = TRUE,
#                          itermax = 10000)
# # save(mrf_seg_shrink, file = "utrecht_mrf_seg.RData")
# load(file = "utrecht_mrf_seg_shrink.RData")
# mrf_seg_shrink$bic %>% plot()
# mrf_seg_shrink$result %>% t() %>% matplot(type = 'l')
# mrf_seg_shrink$result %>% apply(2, function(vec) length(unique(vec)))
# ind <- 7
# polygon["seg_shrinkage"] <- mrf_seg_shrink$result[, ind]
# 
# col_seg_shrink <- div_pal(polygon[["seg_shrinkage"]], by = 0.01)$col
# pal_seg_shrink <- div_pal(polygon[["seg_shrinkage"]], by = 0.01)$pal
# plot(polygon["seg_shrinkage"], col = col_seg_shrink,
#      # main = paste("seg_shrinkage:",
#      #              apply(mrf_seg_shrink$result, 2, function(vec) length(unique(vec)))[ind],
#      #              "params"),
#      main = NA,
#      border = NA)
# legend("right",
#        legend = c("below -0.25", "-0.25 to -0.20", "-0.20 to -0.15", "-0.15 to -0.10",
#                   "-0.10 to -0.05", "-0.05 to 0.00", "0.00 to 0.05", "0.05 to 0.10",
#                   "0.10 to 0.15", "0.15 to 0.20", "0.20 to 0.25", "above 0.25"),
#        fill = pal_seg_shrink(12),
#        bty = "n",
#        cex = 0.58,
#        ncol = 1)
# dev.copy2pdf(file = "utrecht_mrf_seg_shrink.pdf")
# 
# 
# 
# ## Define unions for better borders
# # aaa <- st_union(polygon["seg_shrinkage"], by_feature = TRUE)
# # class(aaa)
# # plot(aaa)
# #
# # aaa <- polygon %>%
# #   st_set_precision(100000) %>%
# #   dplyr::group_by(seg_shrinkage) %>%
# #   dplyr::summarise(geom_fused = st_union(geometry)) %>%
# #   ungroup()
# # plot(aaa)
# # class(aaa)
# # plot(aaa[1])
# # plot(st_geometry(aaa))
# 
# # With Covariance, without shrinkage --------------------------------------
# # mrf_cov_ns <- agraph_prec(gamma = gamma,
# #                           graph = graph,
# #                           prec = prec / 2,
# #                           lambda = 10 ^ seq(-4, -1, length = 100),
# #                           shrinkage = FALSE,
# #                           tol = 1e-8,
# #                           itermax = 50000)
# # save(mrf_cov_ns, file = "utrecht_mrf_cov_ns.RData")
# # load(file = "utrecht_mrf_cov_ns.RData")
# mrf_cov_ns$bic %>% plot()
# mrf_cov_ns$result %>% apply(1, function(vec) length(unique(vec)))
# mrf_cov_ns$result %>% matplot(type = 'S')
# ind <- 1
# polygon["cov_ns"] <- mrf_cov_ns$result[ind, ]
# col_cov_ns <- div_pal(polygon[["cov_ns"]], range = c(-0.25, 0.25), by = 0.01)$col
# pal_cov_ns <- div_pal(polygon[["cov_ns"]], range = c(-0.25, 0.25), by = 0.01)$pal
# plot(polygon["cov_ns"], col = col_cov_ns, main = NA, border = NA)
# legend("right",
#        legend = c("below -0.25", "-0.25 to -0.20", "-0.20 to -0.15", "-0.15 to -0.10",
#                   "-0.10 to -0.05", "-0.05 to 0.00", "0.00 to 0.05", "0.05 to 0.10",
#                   "0.10 to 0.15", "0.15 to 0.20", "0.20 to 0.25", "above 0.25"),
#        fill = pal_cov_ns(12),
#        bty = "n",
#        cex = 0.58,
#        ncol = 1)
# dev.copy2pdf(file = "utrecht_mrf_cov_ns.pdf")
# 
# 
# 
# ## New implementation
# lambda <- 10 ^ seq(-9, 4, length = 50)
# res <- agraph(gamma, graph, #prec = prec,
#               lambda = lambda, itermax = 50000)
# # theta <- res$result[which.min(res$aic),]
# theta <- res$result[23,]
# plot(theta)
# matplot(log10(lambda), res$result, type = "l")
# 
# x.breaks <- c(min(polygon[["gamma"]]),
#               seq(from = -0.25, to = 0.25, by = 0.01),
#               max(polygon[["gamma"]]))
# pal_div <- c(
#   brewer.pal(n = 9, name = "Blues") %>% rev,
#   brewer.pal(n = 9, name = "Oranges")) %>%
#   colorRampPalette
# polygon_fused <- polygon %>%
#   mutate(theta = theta) %>%
#   mutate(group = cut(theta, breaks = sort(unique(theta)),
#                      include.lowest = TRUE)) %>%
#   group_by(group) %>%
#   summarize(theta_fused = mean(theta), do_union = TRUE)
# palette_signal <- pal_div(n = length(x.breaks) - 1)[cut(x = polygon[["gamma"]],
#                                                         breaks = x.breaks, right = FALSE,
#                                                         include.lowest = TRUE, labels = FALSE)]
# palette_seg <- pal_div(n = length(x.breaks) - 1)[cut(x = theta,
#                                                      breaks = x.breaks, right = FALSE,
#                                                      include.lowest = TRUE, labels = FALSE)]
# palette_seg_fused <- pal_div(n = length(x.breaks) - 1)[cut(x = polygon_fused[["theta_fused"]],
#                                                            breaks = x.breaks, right = FALSE,
#                                                            include.lowest = TRUE, labels = FALSE)]
# op <- par(mar = rep(0, 4))
# plot(st_geometry(polygon), col = palette_signal, border = NA)
# dev.copy2pdf(file = "figure/data_utrecht_unseg.pdf", width = 5, height = 5)
# plot(st_geometry(polygon_fused), col = palette_seg_fused, lwd = 2)
# dev.copy2pdf(file = "figure/data_utrecht_seg.pdf", width = 5, height = 5)
# par(op)
