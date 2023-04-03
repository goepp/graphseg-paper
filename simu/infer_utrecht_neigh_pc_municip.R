library(tidyverse); library(gdata); library(RColorBrewer)
library(rgeos); library(rgdal)
library(glmnet); library(igraph); library(sf); library(Matrix)
library(pryr); library(limSolve); library(graphseg)
library(scales); library(latex2exp); library(flsa)
library(kableExtra); library(ggpubr); library(dendextend)

source("utils/sf2nb.R")
source("utils/div_pal.R")
source("utils/infer_functions.R")
width_pdf = 5
height_pdf = 5
pal_div <- (c(brewer.pal(n = 9, name = "Blues") %>% rev,
              brewer.pal(n = 9, name = "Oranges"))) %>% colorRampPalette

utrecht_neigh <- readRDS("simu/synthetic/utrecht_neigh.rds")
graph_utrecht_neigh <- readRDS("simu/synthetic/graph_utrecht_neigh.rds")


# Load data, simulate data ---------------------------------------------
set.seed(87)
baseline <- data.frame("municip_code" = levels(utrecht_neigh$municip_code),
                       "baseline" = rpois(nlevels(utrecht_neigh$municip_code), 10))
utrecht_neigh <- utrecht_neigh %>% mutate(municip_code = fct_drop(municip_code)) %>%
  left_join(baseline, by = "municip_code") %>%
  mutate(true_labels = municip_code %>% as.factor %>% as.numeric)

## Simulate data
nrep <- 100
# sigma_sq <- exp(seq(log(0.05), log(20), length.out = 10))
sigma_sq <- c(seq(0.1, 1.2, 0.1), 2, 3, 5)
names(sigma_sq) <- sigma_sq

set.seed(92)
gamma_array <- replicate(
  nrep, sapply(sigma_sq, function(sigma_sq) utrecht_neigh$baseline + rnorm(nrow(utrecht_neigh), 0, sigma_sq))
)
rms0 <- sapply(seq_along(sigma_sq), function(a)
  sqrt(sum((gamma_array[, a, 1] - utrecht_neigh$baseline) ^ 2 ) / length(gamma))
)


## Perform simulation
lambda <- 10 ^ seq(-3, 2, length = 50)


# Compare (agraph vs flsa) whole regularization paths without rep --------
# Agraph
full_infer_res <- apply(gamma_array[, , 1], 2, agraph, graph = graph_utrecht_neigh,
                        lambda = lambda, weights = NULL, shrinkage = TRUE, tol = 1e-8, itermax = 50000)
rms_path <- lapply(full_infer_res, function(a)
  sqrt(rowSums((a$result - t(replicate(length(lambda), utrecht_neigh$baseline))) ^ 2) / length(gamma)))
dim_path <- lapply(full_infer_res, function(a) a$model_dim)
clust_score_path <- lapply(full_infer_res, function(z) {
  apply(z$result, 1, function(zz) {
    score <- score(as.numeric(as.factor(zz)), true_labels)
    return(ifelse(is.nan(score), 0, score))
  }) %>% unname
})
# Flsa
full_infer_res_flsa <- apply(gamma_array[, , 1], 2, flsa_graph, graph = graph_utrecht_neigh,
                             lambda = lambda)
rms_path_flsa <- lapply(full_infer_res_flsa, function(a)
  sqrt(rowSums((a$result - t(replicate(length(lambda), utrecht_neigh$baseline))) ^ 2) / length(gamma)))
dim_path_flsa <- lapply(full_infer_res_flsa, function(a) a$model_dim)
clust_score_path_flsa <- lapply(full_infer_res_flsa, function(z) {
  apply(z$result, 1, function(zz) {
    score <- score(as.numeric(as.factor(zz)), utrecht_neigh$true_labels)
    return(ifelse(is.nan(score), 0, score))
  }) %>% unname
})

# Plot RMS for a fixed sigma value
criterion_df <- bind_rows(
  bind_cols(list("aic" = lapply(full_infer_res, function(a) a$aic) %>% unlist,
                 "bic" = lapply(full_infer_res, function(a) a$bic) %>% unlist,
                 "gcv" = lapply(full_infer_res, function(a) a$gcv) %>% unlist)) %>%
    dplyr::mutate("lambda_" = rep(lambda, length(sigma_sq)),
                  "Noise_sd" = rep(sigma_sq, each = length(lambda)),
                  "treatment" = "agraph"),
  bind_cols(list("aic" = lapply(full_infer_res, function(a) a$aic) %>% unlist,
                 "bic" = lapply(full_infer_res, function(a) a$bic) %>% unlist,
                 "gcv" = lapply(full_infer_res, function(a) a$gcv) %>% unlist)) %>%
    dplyr::mutate("lambda_" = rep(lambda, length(sigma_sq)),
                  "Noise_sd" = rep(sigma_sq, each = length(lambda)),
                  "treatment" = "flsa"))
path_df <- bind_rows(
  bind_cols(gather(data.frame(dim_path), key = "sigma_sq", value = "model_dim"),
            gather(data.frame(rms_path), key = "sigma_sq", value = "rms"),
            gather(data.frame(clust_score_path), key = "sigma_sq", value = "clust_score")) %>%
    mutate("Noise_sd" = as.numeric(str_remove(sigma_sq...1, "X"))) %>%
    select(-c(sigma_sq...1, sigma_sq...3, sigma_sq...5)) %>%
    mutate("lambda_" = rep(lambda, length(sigma_sq))) %>%
    mutate("treatment" = "agraph"),
  bind_cols(gather(data.frame(dim_path_flsa), key = "sigma_sq", value = "model_dim"),
            gather(data.frame(rms_path_flsa), key = "sigma_sq", value = "rms"),
            gather(data.frame(clust_score_path_flsa), key = "sigma_sq", value = "clust_score")) %>%
    mutate("Noise_sd" = as.numeric(str_remove(sigma_sq...1, "X"))) %>%
    select(-c(sigma_sq...1, sigma_sq...3, sigma_sq...5)) %>%
    mutate("lambda_" = rep(lambda, length(sigma_sq))) %>%
    mutate("treatment" = "flsa")
) %>% inner_join(criterion_df, by = c("Noise_sd", "treatment", "lambda_"))

ggplot() +
  facet_grid(Noise_sd ~., labeller = labeller(Noise_sd = label_both)) +
  geom_line(data = path_df %>% filter(Noise_sd %in% c(0.5, 0.8, 1.0)),
            mapping = aes(model_dim, rms, color = treatment)) +
  geom_point(data = path_df %>% filter(Noise_sd %in% c(0.5, 0.8, 1.0)),
             mapping = aes(model_dim, rms, color = treatment)) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  ylab("Root mean square") + xlab("Model dimension") +
  scale_x_log10() + scale_y_log10() +
  geom_point(data = path_df %>% group_by(treatment, Noise_sd) %>%
               slice(c(which.min(aic), which.min(bic))) %>%
               filter(Noise_sd %in% c(0.5, 0.8, 1.0)) %>% ungroup(),
             mapping = aes(model_dim, rms, colour = treatment),
             size = 4, fill = "white", shape = rep(23:24, 6))
ggsave(file = "simu/figure/utrecht_neigh_pc_municip/utrecht_neigh_pc_municip_rms_model_dim_crit.pdf",
       width = 6, height = 11)

ggplot(path_df %>% filter(Noise_sd %in% c(0.5, 0.8, 1.0)),
       aes(model_dim, rms, color = treatment)) +
  facet_grid(Noise_sd ~., labeller = labeller(Noise_sd = label_both)) +
  geom_line() + geom_point() +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = c(0.15, 0.10)) +
  ylab("Root mean square") + xlab("Model dimension") +
  scale_x_log10() + scale_y_log10()
ggsave(file = "simu/figure/utrecht_neigh_pc_municip/utrecht_neigh_pc_rms_model_dim.pdf",
       width = 6, height = 11)

# Visualizine piecewise constant signal in 1D for aridge and flsa
par(mfrow = c(1, 2))
plot(gamma_array[, ind, 1])
points(full_infer_res_flsa[[ind]]$result[which.min(full_infer_res[[ind]]$aic), ], col = "blue", pch = 3)
points(utrecht_neigh$baseline, col = "red")
plot(gamma_array[, ind, 1])
points(full_infer_res[[ind]]$result[which.min(full_infer_res[[ind]]$aic), ], col = "green", pch = 3)
points(utrecht_neigh$baseline, col = "red")


# Compare with repetition the models selected by AIC, BIC, and G --------
## Inference ---------------------------------------------------------------
# infer_res <- apply(gamma_array, c(2, 3), infer, graph = graph_utrecht_neigh, lambda = lambda,
#                    baseline = utrecht_neigh$baseline,
#                    true_labels = utrecht_neigh$true_labels)
# infer_res_flsa <- apply(gamma_array, c(2, 3), infer_flsa, graph = graph_utrecht_neigh, lambda = lambda,
#                         baseline = utrecht_neigh$baseline, 
#                         true_labels = utrecht_neigh$true_labels)
infer_res <- readRDS(file = "infer_res8.rds")
infer_res_flsa <- readRDS(file = "infer_res_flsa8.rds")


## Post-treatment ----------------------------------------------------------
### Root mean squared error ------------------------------------------------
rms <- apply(infer_res, c(1, 2), function(z) unlist(z[[1]]$rms)) %>%
  "colnames<-"(sigma_sq)
rms_flsa <- apply(infer_res_flsa, c(1, 2), function(z) unlist(z[[1]]$rms)) %>%
  "colnames<-"(sigma_sq)
rms_df <- bind_rows(
  reshape2::melt(rms) %>%
    transmute("criterion" = Var1, "sigma_sq" = as.factor(Var2), "ind_rep" = Var3, "rms" = value) %>%
    mutate("treatment" = "aridge"),
  reshape2::melt(rms_flsa) %>%
    transmute("criterion" = Var1, "sigma_sq" = as.factor(Var2), "ind_rep" = Var3, "rms" = value) %>%
    mutate("treatment" = "flsa")
)

### Dimension ----------------------------------------------------
dim <- apply(infer_res, c(1, 2), function(z) unlist(z[[1]]$dim)) %>%
  "colnames<-"(sigma_sq)
dim_flsa <- apply(infer_res_flsa, c(1, 2), function(z) unlist(z[[1]]$dim)) %>%
  "colnames<-"(sigma_sq)
dim_df <- bind_rows(
  reshape2::melt(dim) %>%
    transmute("criterion" = Var1, "sigma_sq" = as.factor(Var2), "ind_rep" = Var3, "dim" = value) %>%
    mutate("treatment" = "aridge", "criterion" = criterion),
  reshape2::melt(dim_flsa) %>%
    transmute("criterion" = Var1, "sigma_sq" = as.factor(Var2), "ind_rep" = Var3, "dim" = value) %>%
    mutate("treatment" = "flsa"))

### Clustering score ----------------------------------------------------
clust_score <- apply(infer_res, c(1, 2), function(z) unlist(z[[1]]$clust_score)) %>%
  "colnames<-"(sigma_sq)
clust_score_flsa <- apply(infer_res_flsa, c(1, 2), function(z) unlist(z[[1]]$clust_score)) %>%
  "colnames<-"(sigma_sq)
clust_score_df <- bind_rows(
  reshape2::melt(clust_score) %>%
    transmute("criterion" = Var1, "sigma_sq" = as.factor(Var2), "ind_rep" = Var3, "clust_score" = value) %>%
    mutate("treatment" = "aridge", "criterion" = criterion),
  reshape2::melt(clust_score_flsa) %>%
    transmute("criterion" = Var1, "sigma_sq" = as.factor(Var2), "ind_rep" = Var3, "clust_score" = value) %>%
    mutate("treatment" = "flsa"))

## Save results ------------------------------------------------------------
# save(dim_df, rms_df, clust_score_df, infer_res, sigma_sq,
#      file = "simu/simu_results/infer_utrecht_neigh_pc_municip_clust_score.RData")
# load(file = "simu/simu_results/infer_utrecht_neigh_pc_municip.RData")
dim_rms_clust_score_df <- dim_df %>%  inner_join(
  rms_df, by = c("sigma_sq", "criterion", "ind_rep", "treatment")) %>%
  inner_join(clust_score_df, by = c("sigma_sq", "criterion", "ind_rep", "treatment"))
# saveRDS(dim_rms_clust_score_df,
#         file = "simu/simu_results/utrecht_neigh_pc_municip_dim_rms_clust_score_df.rds")
dim_rms_clust_score_df <- readRDS(
  file = "simu/simu_results/utrecht_neigh_pc_municip_dim_rms_clust_score_df.rds")

rms0_df <- tibble(sigma_sq = sigma_sq %>% as.character() %>% as.factor(), 
                  rms0 = rms0)
tidy_df <- dim_rms_clust_score_df %>% left_join(rms0_df, by = "sigma_sq")
saveRDS(tidy_df, file = "simu/simu_results/utrecht_neigh_pc_municip_tidy_df.rds")

## Table results -----------------------------------------------------------
### RMSE and model_dim tables with AIC -------------------------------------
ccc <- dim_rms_df %>%
  filter(criterion == "aic") %>%
  gather(key = "quantity", value = "value", -c(criterion, sigma_sq, ind_rep, treatment)) %>%
  group_by(treatment, sigma_sq, quantity) %>%
  dplyr::summarize(mean = mean(value),
                   sd = sd(value)) %>%
  mutate(mean = signif(mean, 3) %>% as.character,
         sd = signif(sd, 3) %>% as.character) %>%
  mutate(value_err = paste('$', mean, "\\pm", sd, '$'))
min_col_rms_dim <- rbind(ccc %>% select(-c(value_err, sd)) %>%
                           reshape2::acast(sigma_sq ~ quantity + treatment) %>%
                           apply(1, function(a) which.min(a[1:2])),
                         ccc %>% select(-c(value_err, sd)) %>%
                           reshape2::acast(sigma_sq ~ quantity + treatment) %>%
                           apply(1, function(a) which.min(a[3:4]) + 2))
prod_table_rms_dim <- ccc %>% reshape2::acast(sigma_sq ~ quantity + treatment) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  rownames_to_column(var = linebreak("Noise \n standard \n deviation")) %>%
  mutate("no_reg" = paste('$', as.character(signif(rms0, digits = 2)), '$'))
for (ind1 in 1:ncol(min_col_rms_dim)) {
  for (ind2 in 1:nrow(min_col_rms_dim)) {
    prod_table_rms_dim[ind1, min_col_rms_dim[ind2, ind1] + 1] <-
      paste0("\\boldmath{",
             prod_table_rms_dim[ind1, min_col_rms_dim[ind2, ind1] + 1],
             "}")
  }
}
ddd <- prod_table_rms_dim %>%
  slice(seq(1, 15, by = 2)) %>%
  kable("latex", booktabs = TRUE, escape = FALSE,
        col.names = linebreak(c("Noise \n standard \n deviation",
                                rep(c("Adaptive \n Ridge", "FLSA"), 2),
                                "No \n Regula- \n rization"), align = "c")
  ) %>%
  add_header_above(c(" ","Model \n dimension" = 2, "RMS" = 3)) %>%
  kable_styling(latex_options = c("striped", "scale_down"), position = "center")
eee <- prod_table_rms_dim %>% slice(seq(1, 15, by = 2))
saveRDS(eee, file = "simu/table/utrecht_neigh_pc_municip/utrecht_neigh_pc_municip_rms_dim_df.rds")
saveRDS(ddd, file = "simu/table/utrecht_neigh_pc_municip/utrecht_neigh_pc_municip_rms_dim_kable.rds")
ddd %>% cat(file = "simu/table/utrecht_neigh_pc_municip/utrecht_neigh_pc_municip_rms_dim_table.tex")


### RMSE Table -------------------------------------------------------------
aaa <- rms_df %>%
  group_by(treatment, criterion, sigma_sq) %>%
  dplyr::summarize(mean = mean(rms),
                   sd = sd(rms)) %>%
  mutate(mean = sprintf("%.3f", round(mean, digits = 4)),
         sd = sprintf("%.3f", round(sd, digits = 4))) %>%
  mutate(value_err = paste('$', mean, "\\pm", sd, '$'))

min_col <- rbind(aaa %>% select(-c(value_err, sd)) %>%
                   reshape2::acast(sigma_sq ~ criterion + treatment) %>%
                   apply(1, function(a) which.min(a[1:2])),
                 aaa %>% select(-c(value_err, sd)) %>%
                   reshape2::acast(sigma_sq ~ criterion + treatment) %>%
                   apply(1, function(a) which.min(a[3:4]) + 2),
                 aaa %>% select(-c(value_err, sd)) %>%
                   reshape2::acast(sigma_sq ~ criterion + treatment) %>%
                   apply(1, function(a) which.min(a[5:6]) + 4))
prod_table <- aaa %>% reshape2::acast(sigma_sq ~ criterion + treatment) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  rownames_to_column(var = linebreak("Noise \n standard \n deviation")) %>%
  mutate("no_reg" = paste('$', sprintf("%.3f", round(rms0, digits = 4)), '$'))
for (ind1 in 1:ncol(min_col)) {
  for (ind2 in 1:nrow(min_col)) {
    prod_table[ind1, min_col[ind2, ind1] + 1] <- paste0("\\boldmath{",
                                                        prod_table[ind1, min_col[ind2, ind1] + 1],
                                                        "}")
  }
}
bbb <- prod_table %>%
  kable("latex", booktabs = TRUE, escape = FALSE,
        col.names = linebreak(c("Noise \n standard \n deviation",
                                rep(c("Adaptive \n Ridge", "FLSA"), 3),
                                "No \n Regula- \n rization"), align = "c")
  ) %>%
  add_header_above(c(" ","AIC" = 2, "BIC" = 2, "GCV" = "2", " ")) %>%
  kable_styling(latex_options = c("striped", "scale_down"), position = "center")

bbb %>% cat(file = "simu/table/utrecht_neigh/pc/utrecht_neigh_pc_rms_table.tex")

# Plot results
## RMS
# Ardige vs FLSA
rms_df %>% filter(treatment == "aridge" & criterion == "bic" |
                    treatment == "flsa" & criterion == "bic") %>%
  group_by(sigma_sq, treatment) %>%
  dplyr::summarize(n = n()) %>% print(n = Inf)

ggplot(data = rms_df %>% filter(criterion != "gcv"),
       aes(x = as.numeric(as.character(sigma_sq)),
           y = rms, group = interaction(as.numeric(as.character(sigma_sq)), treatment),
           color = treatment)) +
  geom_boxplot() +
  facet_grid(criterion ~.) +
  theme_bw() + ylab("Root Mean Square") + xlab("Noise standard deviation") +
  theme(legend.position = c(0.15, 0.85)) +
  scale_y_log10() + scale_x_log10() +
  stat_summary(fun.y = median, geom = "line", aes(group = treatment))

# previous plot: Aridge only
ggplot(data = rms_df, aes(x = as.numeric(as.character(sigma_sq)), y = rms, group = sigma_sq)) +
  geom_boxplot() +
  facet_grid(criterion ~.) +
  theme_bw() + ylab("Root Mean Square") + xlab("Noise standard deviation") +
  theme(legend.position = c(0.85, 0.17)) +
  scale_y_log10() + scale_x_log10()
ggsave(file = "simu/figure/utr_neigh_pc_rms.pdf", width = 6, height = 5)
## Dimension of the selected model
ggplot(data = dim_df %>% filter(criterion != "gcv"),
       aes(x = sigma_sq, y = dim, color = criterion)) +
  geom_boxplot() +
  theme_bw() + ylab("Effective dimension of the segmentation") + xlab("Noise standard deviation") +
  theme(legend.position = c(0.17, 0.85)) +
  scale_y_log10()
ggsave(file = "simu/figure/utr_neigh_pc_dim.pdf", width = 6, height = 5)
## The model dimension for the GCV increases from 2000 (for sigma_sq = 0.1)
# to 2860 (for sigma = 1)


# Run agraph once and plot results -------
## Using AIC criterion

for (ind_sigma in c(5)) {
  # Perform inference once
  res <- agraph(gamma = gamma_array[, ind_sigma, 1],
                graph = graph_utrecht_neigh,
                lambda = lambda,
                weights = NULL,
                shrinkage = TRUE,
                tol = 1e-8,
                itermax = 50000)

  res_flsa <- flsa_graph(gamma_array[, ind_sigma, 1], graph = graph_utrecht_neigh,
                     lambda = lambda)
  
  utrecht_neigh$theta <- res$result[which.min(res$aic), ]
  utrecht_neigh$theta_flsa <- res_flsa$result[which.min(res$aic), ]

  # res2 <- flsa_graph(gamma_array[, ind_sigma, 1], graph = graph_utrecht_neigh,
  #                    lambda = lambda)
  # utrecht_neigh$theta_flsa <- res2$result[which.min(res$bic), ]

  # Print and save regularization plot
  res_path <- as_tibble(res$result) %>% mutate(lambda = lambda) %>%
    gather(key = "area", value = "value", -lambda)
  res_path_flsa <- as_tibble(res_flsa$result) %>% mutate(lambda = lambda) %>%
    gather(key = "area", value = "value", -lambda)
  ggplot(data = res_path, aes(lambda, value, group = area)) +
    geom_line(size = 0.2) +
    scale_x_log10() +
    xlab(TeX("$\\log_{10}(\\lambda)$")) + ylab("Parameters") +
    theme_minimal() +
    geom_vline(xintercept = lambda[sapply(res[c("aic", "bic", "gcv")], which.min)], col = "red") +
    theme(legend.position = "none")
  # ggsave(file = paste0("simu/figure/utrecht_neigh_pc_municip/utrecht_neigh_pc_municip_reg_path_sigma",
  #                      ind_sigma,".pdf"), width = width_pdf, height = height_pdf)

  ggplot(data = res_path_flsa, aes(lambda, value, group = area)) +
    geom_line(size = 0.2) +
    scale_x_log10() +
    xlab(TeX("$\\log_{10}(\\lambda)$")) + ylab("Parameters") +
    theme_minimal() +
    geom_vline(xintercept = lambda[sapply(res[c("aic", "bic", "gcv")], which.min)], col = "red") +
    theme(legend.position = "none")
  # ggsave(file = paste0("simu/figure/utrecht_neigh_pc_municip/utrecht_neigh_pc_municip_reg_path_flsa_sigma",
  #                      ind_sigma,".pdf"), width = 6.5, height = height_pdf)

  
  utrecht_neigh_fused <- utrecht_neigh %>%
    mutate(group = cut(theta, breaks = sort(unique(theta)),
                       include.lowest = TRUE)) %>%
    group_by(group) %>%
    dplyr::summarize(theta = mean(theta), do_union = TRUE)
  
  utrecht_neigh_fused_flsa <- utrecht_neigh %>%
    mutate(group = cut(theta_flsa, breaks = sort(unique(theta_flsa)),
                       include.lowest = TRUE)) %>%
    group_by(group) %>%
    dplyr::summarize(theta_flsa = mean(theta_flsa), do_union = TRUE)
  
  utrecht_neigh_baseline_fused <- utrecht_neigh %>%
    mutate(group = cut(baseline, breaks = sort(unique(baseline)),
                       include.lowest = TRUE)) %>%
    group_by(group) %>%
    dplyr::summarize(baseline = mean(baseline), do_union = TRUE)
  
  # Plot segmentations without color
  # op <- par(mar = rep(0, 4), mfrow = c(1, 2))
  # plot(st_geometry(utrecht_neigh_fused), main = "", key.pos = NULL)
  # plot(st_geometry(utrecht_neigh), add = TRUE, lwd = 0.2, border = "grey")
  # plot(st_geometry(utrecht_neigh_baseline_fused), main = "", key.pos = NULL)
  # plot(st_geometry(utrecht_neigh), add = TRUE, lwd = 0.2, border = "grey")
  # par(op)
  # dev.copy2pdf(file = paste0("simu/figure/utrecht_neigh/pc/utr_neigh_pc_baseline_vs_seg_sigma",
  #                            ind_sigma,".pdf"), width = 7, height = 5)

  # Define palette of initial signal and baseline signal
  x.breaks <- seq(from = min(gamma_array[, ind_sigma, 1]),
                  to = max(gamma_array[, ind_sigma, 1]),
                  by = 0.01)
  palette_signal <- pal_div(n = length(x.breaks) - 1)[cut(x = gamma_array[, ind_sigma, 1],
                                                          breaks = x.breaks, right = FALSE,
                                                          include.lowest = TRUE, labels = FALSE)]
  palette_baseline_fused <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_neigh_baseline_fused$baseline,
                                                                  breaks = x.breaks, right = FALSE,
                                                                  include.lowest = TRUE, labels = FALSE)]
  # Define palette of segmented signal
  palette_seg <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_neigh$theta,
                                                       breaks = x.breaks, right = FALSE,
                                                       include.lowest = TRUE, labels = FALSE)]
  palette_seg_flsa <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_neigh$theta_flsa,
                                                       breaks = x.breaks, right = FALSE,
                                                       include.lowest = TRUE, labels = FALSE)]
  palette_seg_fused <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_neigh_fused$theta,
                                                             breaks = x.breaks, right = FALSE,
                                                             include.lowest = TRUE, labels = FALSE)]
  palette_seg_fused_flsa <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_neigh_fused_flsa$theta_flsa,
                                                             breaks = x.breaks, right = FALSE,
                                                             include.lowest = TRUE, labels = FALSE)]
  # Plot baseline signal
  op <- par(mar = rep(0, 4))
  plot(st_geometry(utrecht_neigh_baseline_fused), main = "",
       lwd = 1.2, border = "black", col = palette_baseline_fused)
  # plot(st_geometry(utrecht_neigh), col = NA, add = TRUE,
  #      main = "", lwd = 0.1, border = alpha("black", 0.4))
  par(op)
  dev.copy2pdf(file = paste0("simu/figure/utrecht_neigh_pc_municip/utrecht_neigh_pc_municip_baseline_sigma",
                             ind_sigma,".pdf"), width = 4, height = 4)

  # Plot unsegmented signal
  op <- par(mar = rep(0, 4))
  plot(st_geometry(utrecht_neigh), col = palette_signal,
       main = "", lwd = 0.1, border = alpha("black", 0.4))
  par(op)
  dev.copy2pdf(file = paste0("simu/figure/utrecht_neigh_pc_municip/utrecht_neigh_pc_municip_unseg_sigma",
                             ind_sigma,".pdf"), width = 4, height = 4)

  # Plot segmented signal AGRAPH
  op <- par(mar = rep(0, 4))
  plot(st_geometry(utrecht_neigh_fused), lwd = 1.2, border = "black",
       col = palette_seg_fused)
  # plot(st_geometry(utrecht_neigh), add = TRUE, lwd = 0.05, border = alpha("black", 0.3))
  par(op)
  dev.copy2pdf(file = paste0("simu/figure/utrecht_neigh_pc_municip/utrecht_neigh_pc_municip_seg_sigma",
                             ind_sigma,".pdf"), width = 4, height = 4)
  
  # Plot segmented signal FLSA
  op <- par(mar = rep(0, 4))
  plot(st_geometry(utrecht_neigh_fused_flsa), lwd = 1.2, border = "black",
       col = palette_seg_fused_flsa)
  # plot(st_geometry(utrecht_neigh), add = TRUE, lwd = 0.05, border = alpha("black", 0.3))
  par(op)
  dev.copy2pdf(file = paste0("simu/figure/utrecht_neigh_pc_municip/utrecht_neigh_pc_municip_seg_flsa_sigma",
                             ind_sigma,".pdf"), width = 4, height = 4)

  # Plot segmentation overlay AGRAPH
  op <- par(mar = rep(0, 4), mfrow = c(1, 1))
  plot(st_geometry(utrecht_neigh), lwd = 0.1, border = "grey")
  plot(st_geometry(utrecht_neigh_fused), main = "", key.pos = NULL,
       border = alpha("red", 1), add = TRUE, lwd = 1)
  plot(st_geometry(utrecht_neigh_baseline_fused), main = "", key.pos = NULL,
       border = alpha("black", 1), add = TRUE)
  par(op)
  dev.copy2pdf(
    file = paste0("simu/figure/utrecht_neigh_pc_municip/utrecht_neigh_pc_municip_seg_overlay_sigma",
                  ind_sigma,".pdf"), width = 4, height = 4)
  
  # Plot segmentation overlay FLSA
  op <- par(mar = rep(0, 4), mfrow = c(1, 1))
  plot(st_geometry(utrecht_neigh), lwd = 0.1, border = "grey")
  plot(st_geometry(utrecht_neigh_fused_flsa), main = "", key.pos = NULL,
       border = alpha("red", 1), add = TRUE, lwd = 1)
  plot(st_geometry(utrecht_neigh_baseline_fused), main = "", key.pos = NULL,
       border = alpha("black", 1), add = TRUE)
  par(op)
  dev.copy2pdf(
    file = paste0("simu/figure/utrecht_neigh_pc_municip/utrecht_neigh_pc_municip_seg_overlay_flsa_sigma",
                  ind_sigma,".pdf"), width = 4, height = 4)
  
  ## Compute (Adjusted) Rand Index
  clust_agraph <- utrecht_neigh$theta %>% 
    as.factor() %>% 
    as.numeric() # to obtain an integer encoding of the clusters
  clust_flsa <- utrecht_neigh$theta_flsa %>% 
    as.factor() %>% 
    as.numeric() # to obtain an integer encoding of the clusters
  clust_true <- utrecht_neigh$municip_code %>% 
    as.factor() %>% 
    as.numeric()
  
  max(clust_agraph)
  max(clust_flsa)
  max(clust_true)
  library(fossil)
  rand.index(clust_agraph, clust_true)
  rand.index(clust_flsa, clust_true)
  
  adj.rand.index(clust_agraph, clust_true)
  adj.rand.index(clust_flsa, clust_true)
  
}


# Visual inspection: what is the signal-to-noise ratio?
ind = 13
plot(gamma_array[, ind, 1])
points(utrecht_neigh$baseline, col = "red", pch = 3)
plot(gamma_array[, ind, 1])
points(utrecht_neigh$baseline, col = "green", pch = 3)

