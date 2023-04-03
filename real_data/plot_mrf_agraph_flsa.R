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
lambda <- 10 ^ seq(-4, 2, by = 0.05) # length = 121


mrf_seg_agraph <- readRDS(file = "real_data/results/utrecht_mrf_seg_prec.rds")
mrf_seg_flsa <- readRDS(file = "real_data/results/utrecht_mrf_seg_flsa.rds")


ind_agraph <- which.min(mrf_seg_prec$aic)
ind_flsa <- which.min(mrf_seg_flsa$aic)
ind


polygon["seg_agraph"] <- mrf_seg_agraph$result[ind_agraph, ]
polygon["seg_flsa"] <- mrf_seg_flsa$result[ind_flsa, ]


# ## Plot both method results ----------------------------------------------
# col_seg_agraph <- div_pal(polygon[["seg_agraph"]], by = 0.01,
#                           lower = -0.3, upper = 0.4)$col
# pal_seg_agraph <- div_pal(polygon[["seg_agraph"]], by = 0.01,
#                           lower = -0.3, upper = 0.4)$pal
# (legend_breaks <- seq(-0.3, 0.4, length = 11)) # where the legend breaks come from
# 
# legend_breaks_names <- mapply(function(a, b) paste0(a, " to ", b), legend_breaks[-1], legend_breaks[-11])
# 
# plot(polygon["seg_agraph"], col = col_seg_agraph,
#      main = NA,
#      border = NA)
# 
# ## Plot difference between the two methods -------------------------------
# mrf_seg_flsa$model_dim %>% plot
# mrf_seg_agraph$model_dim %>% points
# model_dim <- 50
# 
# ind_flsa <- which(mrf_seg_flsa$model_dim <= model_dim)[1]
# ind_agraph <- which(mrf_seg_agraph$model_dim <= model_dim)[1]
# polygon["diff"] <- mrf_seg_flsa$result[ind_flsa, ] - mrf_seg_agraph$result[ind_agraph, ]
# plot(polygon["diff"])


# Plot agraph & flsa segmentations for 4 dims ------------------------------------
dim_for_paper <- c(700, 457, 158)
for (model_dim in dim_for_paper) {
  ind_flsa <- which(mrf_seg_flsa$model_dim <= model_dim)[1]
  ind_agraph <- which(mrf_seg_agraph$model_dim <= model_dim)[1]
  
  ### Plot agraph
  polygon_fused <- polygon %>%
    mutate(seg = mrf_seg_agraph$result[ind_agraph, ]) %>%
    mutate(group = cut(seg, breaks = sort(unique(seg)),
                       include.lowest = TRUE)) %>%
    group_by(group) %>%
    summarize(seg = mean(seg), do_union = TRUE)
  col_cov_s = div_pal(polygon_fused[["seg"]], by = 0.01, lower = -0.3, upper = 0.4)$col
  op <- par(mar = rep(0, 4))
  plot(polygon_fused["seg"], col = col_cov_s, key.pos = NULL,
       border = "black", # optional: white borders
       lwd = 0.5, # border = NA, # optional: removes borders
       main = NA)
  par(op)  
  dev.copy2pdf(file = paste0("real_data/figure/shrink/mrf_cov_s_",
                             model_dim, ".pdf"),
               width = 7, height = 7)
  
  ### Plot flsa
  polygon_fused <- polygon %>%
    mutate(seg = mrf_seg_flsa$result[ind_flsa, ]) %>%
    mutate(group = cut(seg, breaks = sort(unique(seg)),
                       include.lowest = TRUE)) %>%
    group_by(group) %>%
    summarize(seg = mean(seg), do_union = TRUE)
  col_cov_s = div_pal(polygon_fused[["seg"]], by = 0.01, lower = -0.3, upper = 0.4)$col
  op <- par(mar = rep(0, 4))
  plot(polygon_fused["seg"], col = col_cov_s, key.pos = NULL,
       border = "black", # optional: white borders
       lwd = 0.5, # border = NA, # optional: removes borders
       main = NA)
  par(op)
  
  dev.copy2pdf(file = paste0("real_data/figure/flsa/mrf_cov_flsa_",
                             model_dim, ".pdf"),
               width = 7, height = 7)
}

## Post-treatment: crop pdfs -----------------------------------------------
## Remark: pdfcrop can be downloaded from texlive: 
## Execute `sudo apt-get install texlive-extra-utils` in a terminal in Linux
system('for FILE in real_data/figure/flsa/*.pdf; do pdfcrop "${FILE}" "${FILE}"; done')
system('for FILE in real_data/figure/shrink/*.pdf; do pdfcrop "${FILE}" "${FILE}"; done')


