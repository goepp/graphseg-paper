library(igraph)
library(dplyr)
library(sf)
library(igraph)

g <- set.vertex.attribute(graph_neth_region,
                          "name",
                          value = LETTERS[1:vcount(graph_neth_region)])

neighbors(g, V(g)[1:2])
subcomponent(graph_neth_region, V(graph_neth_region)[1])

subgraph <- function(graph, source, n) {
  vertices_subgraph <- igraph::V(graph)[order(igraph::distances(graph, source))]
  induced <- igraph::induced_subgraph(graph, vertices_subgraph[1:n])
  return(induced)
}

neth_neigh <- readRDS("simu/synthetic/neth_neigh.rds")
graph_neth_neigh <- readRDS("simu/synthetic/graph_neth_neigh.rds")
plot(neth_neigh["province_name"])


set.seed(450)
# Extract a subgraph with 3000 areas -----------------------------------------
number_of_areas <- 3000
# We start with neighborhood BU08281508, which is in a place where there are many
# small neighborhoods.
graph_sub_q <- subgraph(graph_neth_neigh, "BU08281508", number_of_areas)
sub_q <- neth_neigh %>% filter(statcode %in% names(V(graph_sub_q)))

define_q_zones <- function(graph_sub_q, sub_q, number_of_zones) {
  germs <- sample(names(V(graph_sub_q)), number_of_zones, replace = FALSE)
  cluster_index <- distances(graph_sub_q)[, germs] %>% apply(1, function(a) which.min(a))
  cluster_names <- germs[cluster_index]
  sub_q <- sub_q %>% 
    mutate(subgraph_code = cluster_names)
  sub_baseline <- data.frame("subgraph_code" = germs,
                             "baseline" = rpois(number_of_zones, 100))
  return(sub_q %>% left_join(sub_q_10_baseline, by = "subgraph_code"))
}


# Simulate signals with different values of q -----------------------------
set.seed(743)
sub_q_10 <- define_q_zones(graph_sub_q, sub_q, number_of_zones = 10) %>% 
  mutate(signal = baseline + rnorm(., 0, 3))
sub_q_30 <- define_q_zones(graph_sub_q, sub_q, number_of_zones = 30) %>% 
  mutate(signal = baseline + rnorm(., 0, 3))
sub_q_100 <- define_q_zones(graph_sub_q, sub_q, number_of_zones = 100) %>% 
  mutate(signal = baseline + rnorm(., 0, 3))
sub_q_300 <- define_q_zones(graph_sub_q, sub_q, number_of_zones = 300) %>% 
  mutate(signal = baseline + rnorm(., 0, 3))

plot(sub_q_10["signal"])
plot(sub_q_30["signal"])
plot(sub_q_100["signal"])
plot(sub_q_300["signal"])

save(graph_sub_q, sub_q_10, sub_q_30, sub_q_100, sub_q_300,
     file = "simu/synthetic/graph_sub_q.RData")

