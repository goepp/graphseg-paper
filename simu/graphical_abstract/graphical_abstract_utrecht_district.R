library(tidyverse); library(gdata); library(RColorBrewer)
library(rgeos); library(rgdal); library(sf)
library(glmnet); library(igraph); library(sf); library(Matrix)
library(pryr); library(limSolve); library(graphseg)
source("utils/sf2nb.R")
source("utils/div_pal.R")
width_pdf = 6
height_pdf = 6

utrecht_district <- readRDS(file = "simu/synthetic/utrecht_district.rds")
graph_utrecht_district <- readRDS(file = "simu/synthetic/graph_utrecht_district.rds")

# Generate baseline
set.seed(92)
baseline_municip <- data.frame("municip_code" = levels(utrecht_district$municip_code),
                               "baseline_municip" = rpois(nlevels(utrecht_district$municip_code), 10))
utrecht_district <- utrecht_district %>% left_join(baseline_municip, by = "municip_code")


set.seed(92)
sigma_sq <- 0.8
utrecht_district$gamma <- utrecht_district$baseline_municip + rnorm(nrow(utrecht_district), 0, sigma_sq)

## Set palette
x.breaks <- c(min(utrecht_district[["gamma"]]),
              seq(from = 2.6, to = 16.2, by = 0.01),
              max(utrecht_district[["gamma"]]))
pal_div <- (c(brewer.pal(n = 9, name = "Blues") %>% rev,
              brewer.pal(n = 9, name = "Oranges"))) %>% colorRampPalette

# Plot illustrative figure 1
palette_signal <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_district[["gamma"]],
                                                        breaks = x.breaks, right = FALSE,
                                                        include.lowest = TRUE, labels = FALSE)]

# palette_seg <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_district[["theta"]],
#                                                      breaks = x.breaks, right = FALSE,
#                                                      include.lowest = TRUE, labels = FALSE)]
op <- par(mar = rep(0, 4))
# plot(utrecht_district["gamma"], lwd = 0.2, main = "", pal = pal_div, key.pos = NULL)
## I have to choose
plot(utrecht_district["gamma"], lwd = 0.2, main = "", col = palette_signal)
par(op)
dev.off()
dev.copy2pdf(file = "simu/figure/graphical_abstract/graphical_abstract_utrecht_district1.pdf",
             width = width_pdf, height = height_pdf)

### Copy as png -------------------------------------------------------------
library(extrafont)
# font_import(pattern = "lmodern*")
png(filename = "simu/figure/graphical_abstract/graphical_abstract_utrecht_district1.png",
    width = 400, height = 400)
op <- par(mar = rep(0, 4), family = "LM Roman 10")
plot(utrecht_district["gamma"], lwd = 0.4, col = palette_signal, 
     main = "", border = "white")
title("1. Noisy spatial signal", line = -26, font.main = 1, cex = 1.6)
par(op)
dev.off()
system('convert simu/figure/graphical_abstract/graphical_abstract_utrecht_district1.png \\
-trim simu/figure/graphical_abstract/graphical_abstract_utrecht_district1.png')

# Plot Illustrative figure 2
# Possibility one: ugly (not used in the end)
# op <- par(mar = rep(0, 4))
# plot(graph_utrecht_district, vertex.size = 2.3,
#      vertex.frame.color = "white",
#      vertex.color = "black",
#      vertex.label.color = "black",
#      edge.color = "black",
#      vertex.label.cex = 0.001,
#      layout = st_coordinates(st_centroid(utrecht_district)))
# par(op)

# Possibility two:
coord <- st_coordinates(st_centroid(utrecht_district))
adj_municip <- as_adjacency_matrix(graph_utrecht_district, type = "both") %>%
  as("symmetricMatrix") %>% as("dsTMatrix")
edge_list <- data.frame(adj_municip@i + 1, adj_municip@j + 1)
segment_df <- cbind(coord[edge_list[, 1], ], coord[edge_list[, 2], ])
ptmat <- segment_df[, 1:4] %>%
  as.matrix() %>%
  .[2:nrow(.), ]
linesegs <- split(ptmat, 1:nrow(ptmat)) %>%
  lapply(., function(x) {
    x <- matrix(x, nrow = 2, byrow = T)
    x <- st_linestring(x)})
final_sf <- st_sfc(linesegs) %>%
  st_sf('ID' = 1:length(.))

op <- par(mar = rep(0, 4))
plot(st_geometry(utrecht_district), lwd = 0.6, border = "grey")
plot(st_geometry(final_sf), lwd = 0.5, add = TRUE)
plot(st_centroid(utrecht_district), add = TRUE, lwd = 1,
     pch = 21,  bg = palette_signal, col = palette_signal)
par(op)
dev.copy2pdf(file = "simu/figure/graphical_abstract/graphical_abstract_utrecht_district2.pdf", 
             width = width_pdf, height = height_pdf)

### Copy as png -------------------------------------------------------------
library(extrafont)
font_import(pattern = "lmodern*")
png(filename = "simu/figure/graphical_abstract/graphical_abstract_utrecht_district2.png",
    width = 400, height = 400)
op <- par(mar = rep(0, 4), family = "LM Roman 10")
plot(st_geometry(utrecht_district), lwd = 0.6, border = "grey")
plot(st_geometry(final_sf), lwd = 0.5, add = TRUE)
plot(st_centroid(utrecht_district), add = TRUE, lwd = 1,
     pch = 21,  bg = palette_signal, col = palette_signal)
title("2. Adding adjacency graph", line = -26, font.main = 1)
par(op)
dev.off()
system('convert simu/figure/graphical_abstract/graphical_abstract_utrecht_district2.png \\
-trim simu/figure/graphical_abstract/graphical_abstract_utrecht_district2.png')


### START LEGACY CODE
# # Plot illustrative figure 3
# op <- par(mar = rep(0.5, 4))
# plot(graph_utrecht_district, vertex.size = 2.3,
#      vertex.frame.color = palette_signal,
#      vertex.color = palette_signal,
#      vertex.label.color = "black",
#      edge.color = "black",
#      edge.width = 0.5,
#      vertex.label.cex = 0.001,
#      layout = st_coordinates(st_centroid(utrecht_district)))
# par(op)
# dev.copy2pdf(file = "simu/figure/graphical_abstract/graphical_abstract_utrecht_district3.pdf",
#              width = width_pdf, height = height_pdf)
### END LEGACY CODE

# Plot illustrative figure 3 & 4
## Run this section to rerun computations -- START
# res <- agraph(utrecht_district$gamma,
#               graph_utrecht_district,
#               lambda = 10 ^ seq(-5, 4, length = 50),
#               shrinkage = TRUE,
#               tol = 1e-8,
#               itermax = 50000)
# saveRDS(res, file = "simu/graphical_abstract_utrecht_district.rds")
## Run this section to rerun computations -- END
res <- readRDS(file = "simu/graphical_abstract/graphical_abstract_utrecht_district.rds")

theta <- res$result[which.min(res$bic), ]
utrecht_district$theta <- theta
palette_seg <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_district[["theta"]],
                                                     breaks = x.breaks, right = FALSE,
                                                     include.lowest = TRUE, labels = FALSE)]
# Plot the segmented signal on the graph structure
utrecht_district_fused <- utrecht_district %>%
  mutate(group = cut(theta, breaks = sort(unique(theta)),
                     include.lowest = TRUE)) %>%
  group_by(group) %>%
  dplyr::summarize(theta = mean(theta), do_union = TRUE)
palette_fused <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_district_fused[["theta"]],
                                                       breaks = x.breaks, right = FALSE,
                                                       include.lowest = TRUE, labels = FALSE)]

# Plot the segmented signal on the spatial regions
op <- par(mar = rep(0, 4))
plot(st_geometry(utrecht_district_fused), lwd = 1.2, main = "", key.pos = NULL, col = palette_fused)
plot(utrecht_district["geometry"], add = TRUE, main = "", key.pos = NULL, lwd = 0.05, border = "grey30")
plot(st_geometry(utrecht_district_fused), lwd = 1.2, main = "", key.pos = NULL, col = NULL, add = TRUE)
par(op)
dev.copy2pdf(file = "simu/figure/graphical_abstract/graphical_abstract_utrecht_district3.pdf",
             width = width_pdf, height = height_pdf)

### Copy as png -------------------------------------------------------------
library(extrafont)
# font_import(pattern = "lmodern*")
png(filename = "simu/figure/graphical_abstract/graphical_abstract_utrecht_district3.png",
    width = 400, height = 400)
op <- par(mar = rep(0, 4), family = "LM Roman 10")
plot(st_geometry(utrecht_district_fused), lwd = 1.2, main = "", key.pos = NULL, col = palette_fused)
plot(utrecht_district["geometry"], add = TRUE, main = "", key.pos = NULL, lwd = 0.3, border = "grey30")
plot(st_geometry(utrecht_district_fused), lwd = 1.2, main = "", key.pos = NULL, col = NULL, add = TRUE)
title("3. Estimated piece-wise constant signal", line = -26, font.main = 1)
par(op)
dev.off()
system('convert simu/figure/graphical_abstract/graphical_abstract_utrecht_district3.png \\
-trim simu/figure/graphical_abstract/graphical_abstract_utrecht_district3.png')


# Plot the segmented signal and display fused zones
op <- par(mar = rep(0, 4))
plot(st_geometry(utrecht_district_fused), main = "", key.pos = NULL, col = palette_fused)
par(op)
dev.copy2pdf(file = "simu/figure/graphical_abstract/graphical_abstract_utrecht_district4.pdf",
             width = width_pdf, height = height_pdf)

### Copy as png -------------------------------------------------------------
library(extrafont)
# font_import(pattern = "lmodern*")
png(filename = "simu/figure/graphical_abstract/graphical_abstract_utrecht_district4.png",
    width = 400, height = 400)
op <- par(mar = rep(0, 4), family = "LM Roman 10")
plot(st_geometry(utrecht_district_fused), main = "", key.pos = NULL, col = palette_fused)
title("4. Resulting spatial clustering", line = -26, font.main = 1)
par(op)
dev.off()
system('convert simu/figure/graphical_abstract/graphical_abstract_utrecht_district4.png \\
-trim simu/figure/graphical_abstract/graphical_abstract_utrecht_district4.png')


#### Concatenate the 4 PNGs into one --------------------------------------
system("convert +append simu/figure/graphical_abstract/graphical_abstract_utrecht_district1.png \\
simu/figure/graphical_abstract/graphical_abstract_utrecht_district2.png \\
       simu/figure/graphical_abstract/temp1.png")
system("convert +append simu/figure/graphical_abstract/graphical_abstract_utrecht_district3.png \\
simu/figure/graphical_abstract/graphical_abstract_utrecht_district4.png \\
       simu/figure/graphical_abstract/temp2.png")
system("convert -append simu/figure/graphical_abstract/temp1.png \\
       simu/figure/graphical_abstract/temp2.png \\
       simu/figure/graphical_abstract/final_figure.png")


