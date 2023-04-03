library(tidyverse); library(gdata); library(RColorBrewer)
library(rgeos); library(rgdal); library(lattice)
library(glmnet); library(igraph); library(sf)
library(mgcv); library(Matrix)
library(pryr); library(limSolve)
library(RANN); library(car)

source("utils/sf2nb.R")
source("utils/div_pal.R")
source("utils/infer_functions.R")

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

mrf_seg_flsa <- flsa_graph(gamma = gamma, graph = graph, lambda = lambda)
# saveRDS(mrf_seg_flsa, file = "real_data/utrecht_mrf_seg_flsa.rds")
# mrf_seg_flsa <- readRDS(file = "real_data/utrecht_mrf_seg_flsa.rds")

# mrf_seg_flsa$aic %>% plot
# mrf_seg_flsa$bic %>% plot
# mrf_seg_flsa$gcv %>% which.min
# 
# mrf_seg_flsa$model_dim %>% plot
# mrf_seg_flsa$nll %>% plot
# 
# plot(mrf_seg_flsa$model_dim)
# 
# mrf_seg_flsa$gcv %>% plot()
# mrf_seg_flsa$gcv %>% log %>%  plot(x = log10(lambda))
# mrf_seg_flsa$nll %>% plot()
# p <- nrow(polygon)
# plot(log(1 / (p * (1 - mrf_seg_flsa$model_dim / p) ^ 2)))
# 
# mrf_seg_flsa$result  %>% matplot(x = lambda, type = 'l', lwd = 0.5, log = "x")
# mrf_seg_flsa$model_dim
# mrf_seg_flsa$model_dim[which.min(mrf_seg_flsa$gcv)]
# mrf_seg_flsa$model_dim[which.min(mrf_seg_flsa$aic)]
# mrf_seg_flsa$model_dim[which.min(mrf_seg_flsa$bic)]
# # ind <- which.min(mrf_seg_flsa$bic)
# ind <- which.min(mrf_seg_flsa$aic)
# 
# polygon["seg_prec"] <- mrf_seg_flsa$result[ind, ]
# col_seg_prec <- div_pal(polygon[["seg_prec"]], by = 0.01,
#                         lower = -0.3, upper = 0.4)$col
# pal_seg_prec <- div_pal(polygon[["seg_prec"]], by = 0.01,
#                         lower = -0.3, upper = 0.4)$pal
# (legend_breaks <- seq(-0.3, 0.4, length = 11)) # where the legend breaks come from
# legend_breaks_names <- mapply(function(a, b) paste0(a, " to ", b), legend_breaks[-1], legend_breaks[-11])
# plot(polygon["seg_prec"], col = col_seg_prec,
#      # main = paste("seg_prec:",
#      #              apply(mrf_seg_flsa$result, 2, function(vec) length(unique(vec)))[ind],
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
# dev.copy2pdf(file = "utrecht_mrf_seg_flsa.pdf")
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