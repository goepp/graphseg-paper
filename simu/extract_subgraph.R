library(igraph)
library(dplyr)
library(sf)
library(igraph)

load("simu/synthetic/neth_region.RData")
load("simu/synthetic/neth_municip.RData")

plot(neth_district["geometry"])
plot(graph_neth_region)
plot(graph_neth_region)

g <- set.vertex.attribute(graph_neth_region,
                                          "name",
                                          value = LETTERS[1:vcount(graph_neth_region)])

plot(g)

neighbors(g, V(g)[1:2])
subcomponent(graph_neth_region, V(graph_neth_region)[1])

subgraph <- function(graph, source, n) {
  vertices_subgraph <- igraph::V(graph)[order(igraph::distances(graph, source))]
  induced <- igraph::induced_subgraph(graph, vertices_subgraph[1:n])
  return(induced)
}

load("simu/synthetic/neth_neigh.RData")
plot(neth_neigh["province_name"])
neth_neigh %>% filter(province_name == "Limburg")

set.seed(450)
graph_sub_100 <- subgraph(graph_neth_neigh, "BU08810000", 100)
sub_100 <- neth_neigh %>% filter(statcode %in% V(graph_sub_100)$name) %>% 
  mutate(signal = rnorm(nrow(.)))

graph_sub_300 <- subgraph(graph_neth_neigh, "BU08810000", 300)
sub_300 <- neth_neigh %>% filter(statcode %in% V(graph_sub_300)$name) %>% 
  mutate(signal = rnorm(nrow(.)))

graph_sub_1000 <- subgraph(graph_neth_neigh, "BU08810000", 1000)
sub_1000 <- neth_neigh %>% filter(statcode %in% V(graph_sub_1000)$name) %>% 
  mutate(signal = rnorm(nrow(.)))

graph_sub_3000 <- subgraph(graph_neth_neigh, "BU08810000", 3000)
sub_3000 <- neth_neigh %>% filter(statcode %in% V(graph_sub_3000)$name) %>% 
  mutate(signal = rnorm(nrow(.)))
plot(sub_3000)

graph_sub_10000 <- subgraph(graph_neth_neigh, "BU08810000", 10000)
sub_10000 <- neth_neigh %>% filter(statcode %in% V(graph_sub_10000)$name) %>% 
  mutate(signal = rnorm(nrow(.)))

save(graph_sub_100, graph_sub_300, graph_sub_1000, graph_sub_3000, graph_sub_10000,
     sub_100, sub_300, sub_1000, sub_3000, sub_10000, file = "simu/synthetic/graph_sub.RData")

