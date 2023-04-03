library(tidyverse); library(gdata); library(RColorBrewer)
library(rgeos); library(rgdal)
library(glmnet); library(igraph); library(sf); library(Matrix)
library(pryr); library(limSolve); library(graphseg)
library(scales); library(latex2exp); library(flsa)
library(kableExtra); library(dendextend)

source("utils/sf2nb.R")
source("utils/div_pal.R")
source("utils/infer_functions.R")
width_pdf = 5
height_pdf = 5
pal_div <- (c(brewer.pal(n = 9, name = "Blues") %>% rev,
              brewer.pal(n = 9, name = "Oranges"))) %>% colorRampPalette

neth_neigh <- readRDS(file = "simu/synthetic/neth_neigh.rds")
graph_neth_neigh <- readRDS(file = "simu/synthetic/graph_neth_neigh.rds")

set.seed(87)
baseline <- data.frame("municip_code" = levels(neth_neigh$municip_code),
                       "baseline" = rpois(nlevels(neth_neigh$municip_code), 10))
neth_neigh <- neth_neigh %>% left_join(baseline, by = "municip_code") %>% 
  mutate(true_labels = municip_code %>% as.factor %>% as.numeric)

## Simulate data
nrep <- 100
# sigma_sq <- exp(seq(log(0.05), log(20), length.out = 10))
sigma_sq <- c(seq(0.1, 1.2, 0.1), 2, 3, 5)
names(sigma_sq) <- sigma_sq

set.seed(92)
gamma_array <- replicate(
  nrep, sapply(sigma_sq, function(sigma_sq) neth_neigh$baseline + rnorm(nrow(neth_neigh), 0, sigma_sq))
)
rms0 <- sapply(seq_along(sigma_sq), function(a)
  sqrt(sum((gamma_array[, a, 1] - neth_neigh$baseline) ^ 2 ) / length(gamma))
)

## Perform simulation
lambda <- 10 ^ seq(-4, 3, length = 50)

# Compare (agraph vs flsa) whole regularization paths without repetition
full_infer_res <- apply(gamma_array[, , 1], 2, agraph, graph = graph_neth_neigh,
                        lambda = lambda, weights = NULL, shrinkage = TRUE, tol = 1e-8, itermax = 50000)
rms_path <- lapply(full_infer_res, function(a)
  sqrt(rowSums((a$result - t(replicate(length(lambda), neigh$baseline))) ^ 2) / length(gamma)))
dim_path <- lapply(full_infer_res, function(a) a$model_dim)
full_infer_res_flsa <- apply(gamma_array[, , 1], 2, flsa_graph, graph = graph_neth_neigh,
                             lambda = lambda)
rms_path_flsa <- lapply(full_infer_res_flsa, function(a)
  sqrt(rowSums((a$result - t(replicate(length(lambda), neigh$baseline))) ^ 2) / length(gamma)))
dim_path_flsa <- lapply(full_infer_res_flsa, function(a) a$model_dim)

# Plot RMS for a fixed sigma value
rms_path_df <- bind_rows(
  bind_cols(gather(data.frame(dim_path), key = "sigma_sq", value = "model_dim"),
            gather(data.frame(rms_path), key = "sigma_sq", value = "rms")) %>%
    mutate("Noise_sd" = as.numeric(str_remove(sigma_sq, "X"))) %>%
    select(-c(sigma_sq1, sigma_sq)) %>%
    mutate("lambda" = rep(lambda, length(sigma_sq))) %>%
    mutate("treatment" = "aridge"),
  bind_cols(gather(data.frame(dim_path_flsa), key = "sigma_sq", value = "model_dim"),
            gather(data.frame(rms_path_flsa), key = "sigma_sq", value = "rms")) %>%
    mutate("Noise_sd" = as.numeric(str_remove(sigma_sq, "X"))) %>%
    select(-c(sigma_sq1, sigma_sq)) %>%
    mutate("lambda" = rep(lambda, length(sigma_sq))) %>%
    mutate("treatment" = "flsa")
)

ggplot(rms_path_df %>% filter(Noise_sd %in% c(0.5, 0.8, 1.0)),
       aes(lambda, rms, color = treatment)) +
  facet_grid(Noise_sd ~., labeller = labeller(Noise_sd = label_both)) +
  geom_line() + geom_point() +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = c(0.15, 0.10)) +
  ylab("Root mean square") + xlab(TeX("$\\lambda$")) +
  scale_x_log10() + scale_y_log10()
ggsave(file = "simu/figure/neth_neigh/pc/neth_neigh_rms_path_lambda.pdf",
       width = 5, height = 9)
ggplot(rms_path_df %>% filter(Noise_sd %in% c(0.5, 0.8, 1.0)),
       aes(model_dim, rms, color = treatment)) +
  facet_grid(Noise_sd ~., labeller = labeller(Noise_sd = label_both)) +
  geom_line() + geom_point() +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = c(0.15, 0.10)) +
  ylab("Root mean square") + xlab("Model dimension") +
  scale_x_log10() + scale_y_log10()
ggsave(file = "simu/figure/neth_neigh/pc/neth_neigh_rms_path_model_dim.pdf",
       width = 6, height = 11)

plot(dim_path_flsa[[ind]], rms_path_flsa[[ind]], type = "l", col = "green")
lines(dim_path[[ind]], rms_path[[ind]], type = "l", col = "red")
ind = 5
plot(log10(lambda), rms_path_flsa[[ind]], type = "l", col = "green")
lines(log10(lambda), rms_path[[ind]], type = "l", col = "red")
abline(h = sqrt(sum((gamma_array[, ind, 1] - neigh$baseline) ^ 2 ) / length(gamma)))
abline(v = log10(lambda[sapply(full_infer_res_flsa[[ind]]
                               [c("bic", "aic", "gcv")], which.min)]), col = "green")
abline(v = log10(lambda[sapply(full_infer_res[[ind]][c("bic", "aic", "gcv")], which.min)]), col = "red")
legend(1.5, 40, c("aridge", "ridge"), col = c("red", "green"), pch = c(1, 1))

par(mfrow = c(1, 2))
plot(gamma_array[, ind, 1])
# points(full_infer_res_ridge[[ind]]$result[which(lambda == 1), ], col = "blue")
points(full_infer_res_flsa[[ind]]$result[which.min(full_infer_res[[ind]]$aic), ], col = "blue", pch = 3)
points(neigh$baseline, col = "red")
plot(gamma_array[, ind, 1])
points(full_infer_res[[ind]]$result[which.min(full_infer_res[[ind]]$aic), ], col = "green", pch = 3)
points(neigh$baseline, col = "red")

# Compare models selected by aic, bic, gcv with repetition
# infer_res <- apply(gamma_array, c(2, 3), infer, graph = graph_neth_neigh, lambda = lambda,
#                    baseline = neth_neigh$baseline,
#                    true_labels = neth_neigh$true_labels)
# infer_res_flsa <- apply(gamma_array, c(2, 3), infer_flsa, graph = graph_neth_neigh, lambda = lambda,
#                         baseline = neth_neigh$baseline,
#                         true_labels = neth_neigh$true_labels)

infer_res <- readRDS(file = "infer_res6.rds")
# infer_res_flsa <- readRDS(file = "infer_res_flsa6.rds")
# infer_res <- readRDS(file = "infer_res9_1_50.rds")
infer_res_flsa <- readRDS(file = "infer_res_flsa9_1_50.rds")

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
#      file = "simu/simu_results/infer_neth_neigh_pc_municip_dim_rms_clust_score.RData")
# load(file = "simu/simu_results/infer_neth_neigh_pc_municip.RData")
dim_rms_clust_score_df <- dim_df %>%  inner_join(
  rms_df, by = c("sigma_sq", "criterion", "ind_rep", "treatment")) %>%
  inner_join(clust_score_df, by = c("sigma_sq", "criterion", "ind_rep", "treatment"))
saveRDS(dim_rms_clust_score_df,
        file = "simu/simu_results/neth_neigh_pc_municip_dim_rms_clust_score_df.rds")
dim_rms_clust_score_df <- readRDS(
  file = "simu/simu_results/neth_neigh_pc_municip_dim_rms_clust_score_df.rds")

rms0_df <- tibble(sigma_sq = sigma_sq %>% as.character() %>% as.factor(), 
                  rms0 = rms0)
tidy_df <- dim_rms_clust_score_df %>% left_join(rms0_df, by = "sigma_sq")
saveRDS(tidy_df, file = "simu/simu_results/neth_neigh_pc_municip_tidy_df.rds")

## Table results

# Rms and dim table
ccc <- dim_rms_df %>%
  filter(criterion == "bic") %>%
  gather(key = "quantity", value = "value", -c(criterion, sigma_sq, ind_rep, treatment)) %>%
  group_by(treatment, sigma_sq, quantity) %>%
  dplyr::summarize(mean = mean(value),
                   sd = sd(value)) %>%
  mutate(mean = signif(mean, digits = 3) %>% as.character,
         sd = signif(sd, digits = 3) %>% as.character) %>%
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
  mutate("no_reg" = paste('$', signif(rms0, digits = 2) %>% as.character, '$'))
for (ind1 in 1:ncol(min_col_rms_dim)) {
  for (ind2 in 1:nrow(min_col_rms_dim)) {
    prod_table_rms_dim[ind1, min_col_rms_dim[ind2, ind1] + 1] <-
      paste0("\\boldmath{",
             prod_table_rms_dim[ind1, min_col_rms_dim[ind2, ind1] + 1],
             "}")
  }
}
ddd = prod_table_rms_dim %>%
  slice(seq(1, 15, by = 2)) %>%
  kable("latex", booktabs = TRUE, escape = FALSE,
        col.names = linebreak(c("Noise \n standard \n deviation",
                                rep(c("Adaptive \n Ridge", "FLSA"), 2),
                                "No \n Regula- \n rization"), align = "c")
  ) %>%
  add_header_above(c(" ","Model \n dimension" = 2, "RMS" = 3)) %>%
  kable_styling(latex_options = c("striped", "scale_down"), position = "center")
eee = prod_table_rms_dim %>% slice(seq(1, 15, by = 2))
saveRDS(eee, file = "simu/table/neth_neigh_pc_municip/neth_neigh_pc_municip_rms_dim_df.rds")
saveRDS(ddd, file = "simu/table/neth_neigh_pc_municip/neth_neigh_pc_municip_rms_dim_kable.rds")
ddd %>% cat(file = "simu/table/neth_neigh_pc_municip/neth_neigh_pc_municip_rms_dim_table.tex")

# RMS
aaa = rms_df %>%
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
bbb = prod_table %>%
  kable("latex", booktabs = TRUE, escape = FALSE,
        col.names = linebreak(c("Noise \n standard \n deviation",
                                rep(c("Adaptive \n Ridge", "FLSA"), 3),
                                "No \n Regula- \n rization"), align = "c")
  ) %>%
  add_header_above(c(" ","AIC" = 2, "BIC" = 2, "GCV" = "2", " ")) %>%
  kable_styling(latex_options = c("striped", "scale_down"), position = "center")

bbb %>% cat(file = "simu/table/neth_neigh_pc_municip/neth_neigh_pc_municip_rms_table.tex")

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
ggsave(file = "simu/figure/neth_neigh_pc_rms.pdf", width = 6, height = 5)
## Dimension of the selected model
ggplot(data = dim_df %>% filter(criterion != "gcv"),
       aes(x = sigma_sq, y = dim, color = criterion)) +
  geom_boxplot() +
  theme_bw() + ylab("Effective dimension of the segmentation") + xlab("Noise standard deviation") +
  theme(legend.position = c(0.17, 0.85)) +
  scale_y_log10()
ggsave(file = "simu/figure/neth_neigh_pc_dim.pdf", width = 6, height = 5)
## The model dimension for the GCV increases from 2000 (for sigma_sq = 0.1)
# to 2860 (for sigma = 1)


## Run agraph once and plot results

for (ind_sigma in seq_along(sigma_sq)) {
  # Perform inference once
  res <- agraph(gamma = gamma_array[, ind_sigma, 1],
                graph = graph_neigh,
                lambda = lambda,
                weights = NULL,
                shrinkage = TRUE,
                tol = 1e-8,
                itermax = 50000)
  res2 <- flsa_graph(gamma_array[, ind_sigma, 1], graph_neigh, lambda)
  
  neigh$theta <- res$result[which.min(res$bic), ]
  neigh$theta_flsa <- res2$result[which.min(res$bic), ]
  
  
  # Print and save regularization plot
  res_path <- as_tibble(res$result) %>% mutate(lambda = lambda) %>%
    gather(key = "area", value = "value", -lambda)
  res_path_flsa <- as_tibble(res2) %>% mutate(lambda = lambda) %>%
    gather(key = "area", value = "value", -lambda)
  ggplot(data = res_path, aes(lambda, value, group = area)) +
    geom_line(size = 0.2) +
    scale_x_log10() +
    xlab(TeX("$\\log_{10}(\\lambda)$")) + ylab("Parameters") +
    theme_minimal() +
    geom_vline(xintercept = lambda[sapply(res[c("aic", "bic", "gcv")], which.min)], col = "red") +
    theme(legend.position = "none")
  ggsave(file = paste0("simu/figure/neth_neigh/pc/neth_neigh_pc_reg_path_sigma",
                       ind_sigma,".pdf"), width = width_pdf, height = height_pdf)
  
  ggplot(data = res_path_flsa, aes(lambda, value, group = area)) +
    geom_line(size = 0.2) +
    scale_x_log10() +
    xlab(TeX("$\\log_{10}(\\lambda)$")) + ylab("Parameters") +
    theme_minimal() +
    geom_vline(xintercept = lambda[sapply(res[c("aic", "bic", "gcv")], which.min)], col = "red") +
    theme(legend.position = "none")
  ggsave(file = paste0("simu/figure/neth_neigh/pc/neth_neigh_pc_reg_path_flsa_sigma",
                       ind_sigma,".pdf"), width = 6.5, height = height_pdf)
  
  # matplot(log10(lambda), res$result, type = "l", lwd = 0.3, lty = 1,
  #         col = "black", ylab = "Parameter values", xlab = TeX("$\\log_{10}(\\lambda)$"))
  # abline(v = log10(lambda)[sapply(res[c("aic", "bic", "gcv")], which.min)],
  #        col = "red")
  
  neigh_fused <- neigh %>%
    mutate(group = cut(theta, breaks = sort(unique(theta)),
                       include.lowest = TRUE)) %>%
    group_by(group) %>%
    dplyr::summarize(theta = mean(theta), do_union = TRUE)
  neigh_baseline_fused <- neigh %>%
    mutate(group = cut(baseline, breaks = sort(unique(baseline)),
                       include.lowest = TRUE)) %>%
    group_by(group) %>%
    dplyr::summarize(baseline = mean(baseline), do_union = TRUE)
  
  # Plot segmentations without color
  # op <- par(mar = rep(0, 4), mfrow = c(1, 2))
  # plot(st_geometry(neigh_fused), main = "", key.pos = NULL)
  # plot(st_geometry(neigh), add = TRUE, lwd = 0.2, border = "grey")
  # plot(st_geometry(neigh_baseline_fused), main = "", key.pos = NULL)
  # plot(st_geometry(neigh), add = TRUE, lwd = 0.2, border = "grey")
  # par(op)
  # dev.copy2pdf(file = paste0("simu/figure/neth_neigh/pc/neth_neigh_pc_baseline_vs_seg_sigma",
  #                            ind_sigma,".pdf"), width = 7, height = 5)
  
  # Plot segmentation overlay
  op <- par(mar = rep(0, 4), mfrow = c(1, 1))
  plot(st_geometry(neigh_fused), main = "", key.pos = NULL,
       border = alpha("red", 0.8))
  plot(st_geometry(neigh_baseline_fused), main = "", key.pos = NULL,
       border = alpha("black", 0.8), add = TRUE)
  plot(st_geometry(neigh), add = TRUE, lwd = 0.2, border = "grey")
  par(op)
  dev.copy2pdf(file = paste0("simu/figure/neth_neigh/pc/neth_neigh_pc_baseline_seg_overlay_sigma",
                             ind_sigma,".pdf"), width = width_pdf, height = height_pdf)
  
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
  palette_seg_fused <- pal_div(n = length(x.breaks) - 1)[cut(x = utrecht_neigh_fused$theta,
                                                             breaks = x.breaks, right = FALSE,
                                                             include.lowest = TRUE, labels = FALSE)]
  
  # Plot unsegmented signal
  utrecht_neigh$gamma <- gamma_array[, ind_sigma, 1]
  op <- par(mar = rep(0, 4))
  plot(utrecht_neigh["gamma"], col = palette_signal,
       main = "", lwd = 0.2, border = alpha("black", 0.6))
  par(op)
  dev.copy2pdf(file = paste0("simu/figure/neth_neigh/pc/neth_neigh_pc_unseg_sigma",
                             ind_sigma,".pdf"), width = width_pdf, height = height_pdf)
  
  # Plot baseline signal
  op <- par(mar = rep(0, 4))
  plot(st_geometry(utrecht_neigh_baseline_fused), main = "",
       lwd = 1.2, border = "black", col = palette_baseline_fused)
  plot(st_geometry(utrecht_neigh), col = NA, add = TRUE,
       main = "", lwd = 0.2, border = alpha("black", 0.6))
  par(op)
  dev.copy2pdf(file = paste0("simu/figure/neth_neigh/pc/neth_neigh_pc_baseline_sigma",
                             ind_sigma,".pdf"), width = width_pdf, height = height_pdf)
  
  # Plot segmented signal
  op <- par(mar = rep(0, 4))
  plot(st_geometry(utrecht_neigh_fused), lwd = 1.2, border = "black",
       col = palette_seg_fused)
  plot(st_geometry(utrecht_neigh), add = TRUE, lwd = 0.2, border = alpha("black", 0.6))
  par(op)
  dev.copy2pdf(file = paste0("simu/figure/neth_neigh/pc/neth_neigh_pc_seg_sigma",
                             ind_sigma,".pdf"), width = width_pdf, height = height_pdf)
}


# Visual inspection: what is the signal-to-noise ratio?
plot(utrecht_neigh$baseline)
points(gamma_array[,5,1], col = "red", pch = 3)
plot(utrecht_neigh$baseline)
points(gamma_array[,8,1], col = "green", pch = 3)

