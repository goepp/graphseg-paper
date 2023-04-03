library(tidyverse); library(gdata); library(RColorBrewer)
library(rgeos); library(rgdal); library(lattice)
library(glmnet); library(igraph); library(sf)
library(mgcv); library(Matrix)
library(pryr); library(limSolve); library(graphseg)
library(graphseg)

source("utils/sf2nb.R")
source("utils/div_pal.R")
width_pdf = 6
height_pdf = 6

# Load packages
library(dplyr)
library(sf)

year <- 2012:2019
region <- "buurt"
region <- "wijk"
region <- "gemeente"
region <- "ggdregio"
baseurl <- "https://geodata.nationaalgeoregister.nl/cbsgebiedsindelingen/wfs?"

for (i in seq_along(year)) {
  region_year <- paste(region, year[i], sep = "_")
  wfsrequest <- paste0("request=GetFeature&outputFormat=json&typeName=cbs_", region_year, "_gegeneraliseerd")
  url <- paste0(baseurl, wfsrequest)
  tmp <- st_read(dsn = url) %>%
    select(statcode, statnaam)
  saveRDS(
    object = tmp,
    file = paste0("simu/cbs_", region, "/cbs_", region_year, ".rds")
    )
}

