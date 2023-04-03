library(tidyverse); library(gdata); library(RColorBrewer)
library(rgeos); library(rgdal)
library(glmnet); library(igraph); library(sf); library(Matrix)
library(pryr); library(limSolve); library(graphseg)
library(scales); library(latex2exp); library(flsa)
library(kableExtra)

source("utils/sf2nb.R")
source("utils/div_pal.R")
width_pdf = 5
height_pdf = 5

## load("synthetic/neth_region.RData") DO NOT USE: regions are not provinces, 
## and we work with provinces

# Parameters -------
op <- par(mar = rep(0, 4))

## Dataset1: neth_municip_pc_province -------
load("simu/synthetic/neth_municip.RData")
load("simu/synthetic/neth_neigh.RData")
neth_region <- neth_neigh %>% 
  group_by(province_code) %>% 
  dplyr::summarize(do_union = TRUE)
plot(st_geometry(neth_municip), border = "grey")
plot(st_geometry(neth_region), add = TRUE, border = "black", lwd = 1)
dev.copy2pdf(file = "simu/figure/dataset1_zones.pdf", width = width_pdf, height = height_pdf)
neth_region$province_code %>% unique %>% length
neth_municip %>% nrow

## Dataset2: utrecht_district_pc_municip -------
load("simu/synthetic/utrecht_district.RData")
utrecht_municip <- utrecht_district %>%
  group_by(municip_code) %>%
  dplyr::summarize(do_union = TRUE)
plot(st_geometry(utrecht_district), border = "grey")
plot(st_geometry(utrecht_municip), add = TRUE, border = "black", lwd = 1)
dev.copy2pdf(file = "simu/figure/dataset2_zones.pdf", width = width_pdf, height = height_pdf)
utrecht_district$statcode %>% unique %>% length
utrecht_municip$municip_code %>% unique %>% length

## Dataset3: utrecht_neigh_pc_province -------
load("simu/synthetic/utrecht_neigh.RData")
utrecht_province <- utrecht_neigh %>%
  group_by(province_code) %>%
  dplyr::summarize(do_union = TRUE)
plot(st_geometry(utrecht_neigh), border = "grey", lwd = 0.3)
plot(st_geometry(utrecht_province), add = TRUE, border = "black", lwd = 1)
dev.copy2pdf(file = "simu/figure/dataset3_zones.pdf", width = width_pdf, height = height_pdf)
utrecht_neigh$statcode %>% unique %>% length
utrecht_province$province_code %>% unique %>% length

## Dataset 4: utrecht_neigh_pc_municip -------
utrecht_municip <- utrecht_neigh %>% 
  group_by(municip_name) %>% 
  summarize(do_union = TRUE)
plot(st_geometry(utrecht_neigh), border = "grey", lwd = 0.3)
plot(st_geometry(utrecht_municip), add = TRUE, border = "black", lwd = 1)
dev.copy2pdf(file = "simu/figure/dataset4_zones.pdf", width = width_pdf, height = height_pdf)
utrecht_neigh$statcode %>% unique %>% length 
utrecht_neigh$municip_code %>% unique %>% length

## Dataset 5: utrecht_neigh_pc_district -------
utrecht_district <- utrecht_neigh %>% 
  group_by(district_code) %>% 
  dplyr::summarize(do_union = TRUE)
plot(st_geometry(utrecht_neigh), border = "grey", lwd = 0.3)
plot(st_geometry(utrecht_district), add = TRUE, border = "black", lwd = 1)
dev.copy2pdf(file = "simu/figure/dataset5_zones.pdf", width = width_pdf, height = height_pdf)
utrecht_neigh$statcode %>% unique %>% length
utrecht_district$district_code %>% unique %>% length

## Dataset 6: neth_neigh_pc_municip -------
load("simu/synthetic/neth_neigh.RData")
neth_district <- neth_neigh %>% 
  group_by(municip_code) %>% 
  dplyr::summarize(do_union = TRUE)
plot(st_geometry(neth_neigh), border = "grey", lwd = 0.3)
plot(st_geometry(neth_district), add = TRUE, border = "black", lwd = 1)
dev.copy2pdf(file = "simu/figure/dataset6_zones.pdf", width = width_pdf, height = height_pdf)
neth_neigh$statcode %>% unique %>% length
neth_district$municip_code %>% unique %>% length

# Parameters -------
par(op)