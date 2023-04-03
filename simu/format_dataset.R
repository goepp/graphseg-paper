library(tidyverse); library(gdata); library(RColorBrewer)
library(rgeos); library(rgdal); library(sf)
library(glmnet); library(igraph); library(sf); library(Matrix)
library(pryr); library(limSolve); library(graphseg)
library(stringi); library(xlsx)

source("utils/sf2nb.R")
source("utils/div_pal.R")
width_pdf = 6
height_pdf = 6

remove_str_as_factor <- function(sf) {
  sf %>%
    as_tibble() %>% # transform to tibble
    st_as_sf() %>% # transform to sf-tibble
    mutate_if(is.factor, as.character) # revert stringsAsFactors
}
prevalence <- read.delim(file = "simu/prev_overgewicht_utrecht.txt")
neth_neigh <- readRDS("simu/cbs_buurt/cbs_buurt_2016.rds") %>% remove_str_as_factor()
neth_district <- readRDS("simu/cbs_wijk/cbs_wijk_2016.rds") %>% remove_str_as_factor()
neth_municip <- readRDS("simu/cbs_gemeente/cbs_gemeente_2016.rds") %>% remove_str_as_factor()
neth_region <- readRDS("simu/cbs_ggdregio/cbs_ggdregio_2016.rds") %>% remove_str_as_factor()


# Plot the different areas ----------------------------------------------------------------------------------------
op <- par(mar = rep(0, 4))
plot(st_geometry(neth_neigh), border = "blue")
plot(st_geometry(neth_district), border = "green", add = TRUE)
plot(st_geometry(neth_municip), add = TRUE)
plot(st_geometry(neth_region), add = TRUE, border = "red")
par(op)

## Number of areas for each dataset
sapply(list(neth_region, neth_municip, neth_district, neth_neigh), nrow)

## Beware: statnaam has duplicates, while statcode does not
## -> use statcode to uniquely identify a row
neth_district %>% count(statcode) %>% filter(n > 1)
neth_district %>% count(statnaam) %>% filter(n > 1)
### same for neth_neigh

## Define province codes for municip

# The file "gemeenten-alfabetisch-2016.xls" comes from the following URL:
# "https://www.cbs.nl/-/media/imported/onze-diensten/methoden/classificaties/documents/2015/39/gemeenten-alfabetisch-2016.xls?la=nl-nl"
# You have to modify the "Gemeentenaam" entry from "Súdwest-Fryslân" to "Sudwest-Fryslan",
# then load it using the following code:
corresp <- read.xlsx("simu/synthetic/gemeenten-alfabetisch-2016.xls", sheetIndex = 1, encoding = "UTF-8",
                     stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  transmute(province_name = Provincienaam %>% stri_trans_general("Latin-ASCII"),
            province_code = Provinciecode %>% as.factor,
            municip_name = Gemeentenaam %>% stri_trans_general("Latin-ASCII"),
            municip_code = Gemeentecode %>% as.factor)

neth_municip <- neth_municip %>%
  mutate(municip_name = statnaam %>% stri_trans_general("Latin-ASCII")) %>%
  left_join(corresp, by = "municip_name") %>%
  mutate(province_name = province_name)

## Define district and municip codes (there are no region codes)
neth_neigh <- neth_neigh %>%
  mutate(district_code = statcode %>% str_remove("BU") %>%
           substr(1, 6) %>% as.factor,
         municip_code = statcode %>% str_remove("BU") %>%
         substr(1, 4) %>% as.factor) %>%
  left_join(corresp, by = "municip_code")
neth_district <- neth_district %>%
  mutate(municip_code = statcode %>% str_remove("WK") %>%
           substr(1, 4) %>% as.factor) %>%
  left_join(corresp, by = "municip_code")

# Define utrecht region
utrecht_neigh <- neth_neigh %>% filter(statcode %in% prevalence$bu_code) %>% 
  mutate_if(is.factor, fct_drop)
utrecht_district <- neth_district %>%
  filter(str_remove(statcode, "WK") %in% (prevalence$bu_code %>% str_remove("BU") %>%
                                                     substr(1,6) %>% unique)) %>%
  mutate(municip_code = statcode %>% as.character() %>% str_remove("WK") %>%
           substr(1, 4) %>% as.factor) %>% 
  mutate(province_code = province_code %>% fct_drop())

# Create adj graph
sf2graph <-
sf2graph <- function(sf) {
  sf %>%
    sf2nb(queen = FALSE, sparse = FALSE) %>%
    'dimnames<-'(list(sf$statcode, sf$statcode)) %>%
    igraph::graph_from_adjacency_matrix(., mode = "undirected")
}
## For municipality subdivision
graph_neth_municip <- neth_municip %>% sf2graph()
## For district subdivision
graph_neth_district <- neth_district %>% sf2graph()
## For neighborhood subdivision
graph_neth_neigh <- neth_neigh %>% sf2graph()
## For region subdivision
graph_neth_region <- neth_region %>% sf2graph()
## For utrecht neigh and district
graph_utrecht_neigh <- utrecht_neigh %>% sf2graph()
graph_utrecht_district <- utrecht_district %>% sf2graph()

# Save data
saveRDS(neth_region, file = "simu/synthetic/neth_region.rds")
saveRDS(graph_neth_region, file = "simu/synthetic/graph_neth_region.rds")

saveRDS(neth_municip, file = "simu/synthetic/neth_municip.rds")
saveRDS(graph_neth_municip, file = "simu/synthetic/graph_neth_municip.rds")

saveRDS(neth_district, file = "simu/synthetic/neth_district.rds")
saveRDS(graph_neth_district, file = "simu/synthetic/graph_neth_district.rds")

saveRDS(neth_neigh, file = "simu/synthetic/neth_neigh.rds")
saveRDS(graph_neth_neigh, file = "simu/synthetic/graph_neth_neigh.rds")

saveRDS(utrecht_neigh, file = "simu/synthetic/utrecht_neigh.rds")
saveRDS(graph_utrecht_neigh, file = "simu/synthetic/graph_utrecht_neigh.rds")

saveRDS(utrecht_district, file = "simu/synthetic/utrecht_district.rds")
saveRDS(graph_utrecht_district, file = "simu/synthetic/graph_utrecht_district.rds")
