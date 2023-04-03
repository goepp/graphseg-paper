library(tidyverse); library(gdata); library(RColorBrewer)
library(rgeos); library(rgdal); library(lattice)
library(glmnet); library(igraph); library(sf)
library(mgcv); library(Matrix)
library(pryr); library(limSolve); library(graphseg)
library(RANN); library(graphseg)
source("utils/sf2nb.R")
source("utils/div_pal.R")

width_pdf <- 6
height_pdf <- 6

# Load raw data -----------------------------------------------------------
mrf <- read.delim(file = "real_data/raw_data/mrfterm_overweight_utrecht.txt")
gamma <- mrf$mrf
sigma_sq <- mrf$se.mrf
polygon <- st_read(dsn = "real_data/raw_data/map_netherland.geojson", quiet = TRUE) %>%
  filter(statcode %in% mrf$bu_code) %>%
  mutate(gamma = gamma)


# Load segmented utrecht mrf signal ---------------------------------------
mrf_seg_prec <- readRDS(file = "real_data/utrecht_mrf_seg_prec.rds")

### Palette and colors
breaks = div_pal(polygon[["gamma"]], by = 0.01, lower = -0.3, upper = 0.4)$breaks
plot(polygon["gamma"], border = NA, col = col_unseg, main = NA)
(legend_breaks <- seq(-0.3, 0.4, length = 11)) # where the legend breaks come from
legend_breaks_names <- mapply(function(a, b) paste0(a, " to ", b), legend_breaks[-1], legend_breaks[-11])

ind <- 70
mrf_seg_prec$model_dim[ind]

ind_for_paper <- c(52, 57, 61, 73)
mrf_seg_prec$model_dim[ind_for_paper] %>% 
  round(digits = 2)

mrf_seg_prec$gcv %>% plot
mrf_seg_prec$aic %>% plot
mrf_seg_prec$bic %>% plot

# desired_model_dims

for (ind in ind_for_paper) {
  polygon_fused <- polygon %>%
    mutate(seg = mrf_seg_prec$result[ind, ]) %>%
    mutate(group = cut(seg, breaks = sort(unique(seg)),
                       include.lowest = TRUE)) %>%
    group_by(group) %>%
    summarize(seg = mean(seg), do_union = TRUE)
  col_cov_s = div_pal(polygon_fused[["seg"]], by = 0.01, lower = -0.3, upper = 0.4)$col
  plot(polygon_fused["seg"], col = col_cov_s, key.pos = NULL,
       border = "black", # optional: white borders
       lwd = 0.5, # border = NA, # optional: removes borders
       main = NA)
  dev.copy2pdf(file = paste0("real_data/figure/shrink/mrf_cov_s_",
                             round(mrf_seg_prec$model_dim[ind]),
                             ".pdf"),
               width = 7, height = 7)
  
  # polygon[paste0("seg_", ind)] <- mrf_seg_prec$result[ind, ]
  # col_seg_prec <- pal_div(n = length(x.breaks) - 1)[cut(x = polygon[[paste0("seg_", ind)]],
  #                                                       breaks = x.breaks, right = FALSE,
  #                                                       include.lowest = TRUE, labels = FALSE)]
  # plot(polygon[paste0("seg_", ind)], col = col_seg_prec,
  #      main = NA,
  #      border = NA,
  #      key.pos = NULL)
  # dev.copy2pdf(file = paste0("real_data/figure/shrink/mrf_cov_s_",
  #                            apply(mrf_seg_prec$result, 2, function(vec) length(unique(vec)))[ind],
  #                            ".pdf"),
  #              width = 7, height = 7)
}


