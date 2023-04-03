library(tidyverse); library(gdata); library(RColorBrewer)
library(rgeos); library(rgdal); library(sf)
library(glmnet); library(igraph); library(sf); library(Matrix)
library(pryr); library(limSolve); library(graphseg)

source("utils/sf2nb.R")
source("utils/div_pal.R")
width_pdf = 6
height_pdf = 6

## load data ---------------------------------------------------------------------
utrecht_neigh <- readRDS(file = "simu/synthetic/utrecht_neigh.rds")
utrecht_district <- readRDS(file = "simu/synthetic/utrecht_district.rds")
neth_municip <- readRDS(file = "simu/synthetic/neth_municip.rds")
neth_neigh <- readRDS(file = "simu/synthetic/neth_neigh.rds")

## utrecht_neigh_pc_municip ---------------------------------------------------------------------
set.seed(87)
baseline_pc <- data.frame("municip_code" = unique(utrecht_neigh$municip_code),
                          "baseline_municip" = rpois(length(unique(utrecht_neigh$municip_code)), 10))
utrecht_neigh <- utrecht_neigh %>% left_join(baseline_pc, by = "municip_code")
sigma_sq <- 1.2
utrecht_neigh$gamma <- utrecht_neigh$baseline_municip + rnorm(nrow(utrecht_neigh), 0, sigma_sq)

## Plot
x.breaks <- c(min(utrecht_neigh[["gamma"]]),
              seq(from = 2.6, to = 16.2, by = 0.01),
              max(utrecht_neigh[["gamma"]]))
pal_div <- (c(brewer.pal(n = 9, name = "Blues") %>% rev,
              brewer.pal(n = 9, name = "Oranges"))) %>% colorRampPalette
palette_signal <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_neigh[["gamma"]],
                                                        breaks = x.breaks, right = FALSE,
                                                        include.lowest = TRUE, labels = FALSE)]
utrecht_neigh_fused <- utrecht_neigh %>%
  group_by(municip_code) %>%
  dplyr::summarize(baseline_municip = mean(baseline_municip), do_union = TRUE)
op <- par(mar = rep(0, 4))
plot(st_geometry(utrecht_neigh), main = "", col = palette_signal, 
     key.pos = NULL, lwd = 0.1, border = NA)
plot(st_geometry(utrecht_neigh_fused), main = "", key.pos = NULL, add = TRUE, lwd = 2)
par(op)
dev.copy2pdf(file = "simu/figure/dataset1.pdf", width = width_pdf, height = height_pdf)

## utrecht_district_pc_municip ---------------------------------------------
set.seed(87)
baseline_pc <- data.frame("municip_code" = unique(utrecht_district$municip_code),
                          "baseline_municip" = rpois(length(unique(utrecht_district$municip_code)), 10))
utrecht_district <- utrecht_district %>% left_join(baseline_pc, by = "municip_code")
sigma_sq <- 1.2
utrecht_district$gamma <- utrecht_district$baseline_municip + rnorm(nrow(utrecht_district), 0, sigma_sq)

## Plot
x.breaks <- c(min(utrecht_district[["gamma"]]),
              seq(from = 2.6, to = 16.2, by = 0.01),
              max(utrecht_district[["gamma"]]))
pal_div <- (c(brewer.pal(n = 9, name = "Blues") %>% rev,
              brewer.pal(n = 9, name = "Oranges"))) %>% colorRampPalette
palette_signal <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_district[["gamma"]],
                                                        breaks = x.breaks, right = FALSE,
                                                        include.lowest = TRUE, labels = FALSE)]
utrecht_district_fused <- utrecht_district %>%
  group_by(municip_code) %>%
  dplyr::summarize(baseline_municip = mean(baseline_municip), do_union = TRUE)
op <- par(mar = rep(0, 4))
plot(st_geometry(utrecht_district), main = "", key.pos = NULL, 
     col = palette_signal, lwd = 0.1, border = NA)
plot(st_geometry(utrecht_district_fused), main = "", key.pos = NULL, add = TRUE, lwd = 2)
par(op)
dev.copy2pdf(file = "simu/figure/dataset2.pdf", width = width_pdf, height = height_pdf)

# utrecht_neigh_pc_province -----------------------------------------------
set.seed(87)
baseline_pc <- data.frame("province_code" = unique(utrecht_neigh$province_code),
                          "baseline_province" = rpois(length(unique(utrecht_neigh$province_code)), 10))
utrecht_neigh <- utrecht_neigh %>% left_join(baseline_pc, by = "province_code")
sigma_sq <- 1.2
utrecht_neigh$gamma <- utrecht_neigh$baseline_province + rnorm(nrow(utrecht_neigh), 0, sigma_sq)

## Plot
x.breaks <- c(min(utrecht_neigh[["gamma"]]),
              seq(from = 2.6, to = 16.2, by = 0.01),
              max(utrecht_neigh[["gamma"]]))
pal_div <- (c(brewer.pal(n = 9, name = "Blues") %>% rev,
              brewer.pal(n = 9, name = "Oranges"))) %>% colorRampPalette
palette_signal <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_neigh[["gamma"]],
                                                        breaks = x.breaks, right = FALSE,
                                                        include.lowest = TRUE, labels = FALSE)]
utrecht_neigh_fused <- utrecht_neigh %>%
  group_by(province_code) %>%
  dplyr::summarize(baseline_province = mean(baseline_province), do_union = TRUE)
op <- par(mar = rep(0, 4))
plot(st_geometry(utrecht_neigh), main = "", key.pos = NULL, 
     col = palette_signal, lwd = 0.1, border = NA)
plot(st_geometry(utrecht_neigh_fused), main = "", key.pos = NULL, add = TRUE, lwd = 2)
par(op)
dev.copy2pdf(file = "simu/figure/dataset3.pdf", width = width_pdf, height = height_pdf)

# utrecht_neigh_pc_district -----------------------------------------------
set.seed(87)
baseline_pc <- data.frame("district_code" = unique(utrecht_neigh$district_code),
                          "baseline_district" = rpois(length(unique(utrecht_neigh$district_code)), 10))
utrecht_neigh <- utrecht_neigh %>% left_join(baseline_pc, by = "district_code")
sigma_sq <- 1.2
utrecht_neigh$gamma <- utrecht_neigh$baseline_district + rnorm(nrow(utrecht_neigh), 0, sigma_sq)

## Plot
x.breaks <- c(min(utrecht_neigh[["gamma"]]),
              seq(from = 2.6, to = 16.2, by = 0.01),
              max(utrecht_neigh[["gamma"]]))
pal_div <- (c(brewer.pal(n = 9, name = "Blues") %>% rev,
              brewer.pal(n = 9, name = "Oranges"))) %>% colorRampPalette
palette_signal <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_neigh[["gamma"]],
                                                        breaks = x.breaks, right = FALSE,
                                                        include.lowest = TRUE, labels = FALSE)]
utrecht_neigh_fused <- utrecht_neigh %>%
  group_by(district_code) %>%
  dplyr::summarize(baseline_district = mean(baseline_district), do_union = TRUE)
op <- par(mar = rep(0, 4))
plot(st_geometry(utrecht_neigh), main = "", key.pos = NULL, 
     col = palette_signal, lwd = 0.05, border = NA)
plot(st_geometry(utrecht_neigh_fused), main = "", key.pos = NULL, add = TRUE, lwd = 2)
par(op)
dev.copy2pdf(file = "simu/figure/dataset4.pdf", width = width_pdf, height = height_pdf)

# neth_neigh_pc_municip ---------------------------------------------------
set.seed(87)
baseline_pc <- data.frame("municip_code" = unique(neth_neigh$municip_code),
                          "baseline_municip" = rpois(length(unique(neth_neigh$municip_code)), 10))
neth_neigh <- neth_neigh %>% left_join(baseline_pc, by = "municip_code")
sigma_sq <- 1.2
neth_neigh$gamma <- neth_neigh$baseline_municip + rnorm(nrow(neth_neigh), 0, sigma_sq)

## Plot
x.breaks <- c(min(neth_neigh[["gamma"]]),
              seq(from = 2.6, to = 16.2, by = 0.01),
              max(neth_neigh[["gamma"]]))
pal_div <- (c(brewer.pal(n = 9, name = "Blues") %>% rev,
              brewer.pal(n = 9, name = "Oranges"))) %>% colorRampPalette
palette_signal <- pal_div(n = length(x.breaks) - 1)[cut(x = neth_neigh[["gamma"]],
                                                        breaks = x.breaks, right = FALSE,
                                                        include.lowest = TRUE, labels = FALSE)]
neth_neigh_fused <- neth_neigh %>%
  group_by(municip_code) %>%
  dplyr::summarize(baseline_municip = mean(baseline_municip), do_union = TRUE)
op <- par(mar = rep(0, 4))
plot(st_geometry(neth_neigh), main = "", key.pos = NULL, 
     col = palette_signal, lwd = 0.05, border = NA)
plot(st_geometry(neth_neigh_fused), main = "", key.pos = NULL, add = TRUE, lwd = 2)
par(op)
dev.copy2pdf(file = "simu/figure/dataset5.pdf", width = width_pdf, height = height_pdf)

## neth_municip_pc_province ---------------------------------------------------------------------
set.seed(87)
baseline_pc <- data.frame("province_code" = unique(neth_municip$province_code),
                          "baseline_province" = rpois(length(unique(neth_municip$province_code)), 10))
neth_municip <- neth_municip %>% left_join(baseline_pc, by = "province_code")
sigma_sq <- 1.2
neth_municip$gamma <- neth_municip$baseline_province + rnorm(nrow(neth_municip), 0, sigma_sq)
## Plot
x.breaks <- c(min(neth_municip[["gamma"]]),
              seq(from = 2.6, to = 16.2, by = 0.01),
              max(neth_municip[["gamma"]]))
pal_div <- (c(brewer.pal(n = 9, name = "Blues") %>% rev,
              brewer.pal(n = 9, name = "Oranges"))) %>% colorRampPalette
palette_signal <- pal_div(n = length(x.breaks) - 1)[cut(x = neth_municip[["gamma"]],
                                                        breaks = x.breaks, right = FALSE,
                                                        include.lowest = TRUE, labels = FALSE)]
neth_municip_fused <- neth_municip %>%
  group_by(province_code) %>%
  dplyr::summarize(baseline_province = mean(baseline_province), do_union = TRUE)
op <- par(mar = rep(0, 4))
plot(st_geometry(neth_municip), lwd = 0.1, main = "", key.pos = NULL, 
     col = palette_signal, border = NA)
plot(st_geometry(neth_municip_fused), main = "", key.pos = NULL, add = TRUE, lwd = 2)
par(op)
dev.copy2pdf(file = "simu/figure/dataset6.pdf", width = width_pdf, height = height_pdf)

