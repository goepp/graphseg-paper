library(tidyverse); library(gdata); library(RColorBrewer)
library(rgeos); library(rgdal)
library(glmnet); library(igraph); library(sf); library(Matrix)
library(pryr); library(limSolve); library(graphseg)
library(scales); library(latex2exp); library(flsa)
library(kableExtra); library(parallel)

source("utils/sf2nb.R")
source("utils/div_pal.R")
source("utils/infer_functions.R")
width_pdf = 5
height_pdf = 5
pal_div <- (c(brewer.pal(n = 9, name = "Blues") %>% rev,
              brewer.pal(n = 9, name = "Oranges"))) %>% colorRampPalette

# load spatial data ---------------------------------------------------------------
region <- "utrecht"
area <- "neigh"
zone <- "district"
zone_code <- paste0(zone, "_code")

dataset <- readRDS(file = paste0("simu/synthetic/", region, "_", area, ".rds"))
graph <- readRDS(file = paste0("simu/synthetic/graph_", region, "_", area, ".rds"))


# load simulation results ---------------------------------------------------------------
full_infer_res <- readRDS(file = paste0("simu/simu_results/", region, "_", area, "_", zone, "_full_infer_res.rds"))
full_infer_res_flsa <- readRDS(file = paste0("simu/simu_results/", region, "_", area, "_", zone, "_full_infer_res_flsa.rds"))


# Plot RMS for a fixed sigma value
criterion_df <- bind_rows(
  bind_cols(list("aic" = lapply(full_infer_res, function(a) a$aic) %>% unlist,
                 "bic" = lapply(full_infer_res, function(a) a$bic) %>% unlist,
                 "gcv" = lapply(full_infer_res, function(a) a$gcv) %>% unlist)) %>%
    dplyr::mutate("treatment" = "aridge",
                  "lambda_" = rep(rep(lambda, length(sigma_sq)), nrep),
                  "ind_rep" = rep(1:nrep, each = length(lambda) * length(sigma_sq)),
                  "Noise_sd" = rep(rep(sigma_sq, each = length(lambda)), nrep)),
  bind_cols(list("aic" = lapply(full_infer_res_flsa, function(a) a$aic) %>% unlist,
                 "bic" = lapply(full_infer_res_flsa, function(a) a$bic) %>% unlist,
                 "gcv" = lapply(full_infer_res_flsa, function(a) a$gcv) %>% unlist)) %>%
    dplyr::mutate("treatment" = "flsa",
                  "lambda_" = rep(rep(lambda, length(sigma_sq)), nrep),
                  "ind_rep" = rep(1:nrep, each = length(lambda) * length(sigma_sq)),
                  "Noise_sd" = rep(rep(sigma_sq, each = length(lambda)), nrep)))
path_df <- bind_rows(
  bind_cols(list(
    "rms" = lapply(full_infer_res, function(a)
      sqrt(rowSums((a$result - t(replicate(length(lambda), dataset$baseline))) ^ 2) / length(gamma))) %>% 
      unlist(), 
    "model_dim" = lapply(full_infer_res, function(a) a$model_dim) %>% unlist())) %>% 
    dplyr::mutate("Noise_sd" = rep(rep(sigma_sq, each = length(lambda)), nrep),
                  "lambda_" = rep(rep(lambda, length(sigma_sq)), nrep),
                  "ind_rep" = rep(1:nrep, each = length(lambda) * length(sigma_sq)),
                  "treatment" = "aridge"),
  bind_cols(list(
    "rms" = lapply(full_infer_res_flsa, function(a)
      sqrt(rowSums((a$result - t(replicate(length(lambda), dataset$baseline))) ^ 2) / length(gamma))) %>% 
      unlist(), 
    "model_dim" = lapply(full_infer_res_flsa, function(a) a$model_dim) %>% unlist())) %>% 
    dplyr::mutate("Noise_sd" = rep(rep(sigma_sq, each = length(lambda)), nrep),
                  "lambda_" = rep(rep(lambda, length(sigma_sq)), nrep),
                  "ind_rep" = rep(1:nrep, each = length(lambda) * length(sigma_sq)),
                  "treatment" = "flsa")) 
all_df <- dplyr::inner_join(path_df, criterion_df, 
                            by = c("Noise_sd", "treatment", "lambda_", "ind_rep"))
reps_out_df <- all_df %>% 
  group_by(Noise_sd, treatment, lambda_) %>% 
  dplyr::summarize_at(vars(rms, model_dim), list(mean = mean, sd = sd)) %>% 
  dplyr::mutate(rms_lower = rms_mean - rms_sd, 
                rms_upper = rms_mean + rms_sd,
                model_dim_lower = model_dim_mean - model_dim_sd, 
                model_dim_upper = model_dim_mean + model_dim_sd)
reps_out_df


# plot rms vs model_dim ---------------------------------------------------
ggplot(reps_out_df %>% filter(Noise_sd %in% c(0.5, 0.8, 1.0)),
       aes(model_dim_mean, rms_mean, color = treatment)) +
  facet_grid(Noise_sd ~ ., labeller = labeller(Noise_sd = label_both),
             scales = "free_y") +
  geom_line(alpha = 0.3) +
  geom_linerange(aes(ymin = rms_lower, ymax = rms_upper)) +
  geom_linerange(aes(xmin = model_dim_lower, xmax = model_dim_upper)) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  ylab("Root mean square") + xlab("Model dimension") +
  scale_y_log10()

ggplot() +
  facet_grid(Noise_sd ~ ., labeller = labeller(Noise_sd = label_both)) +
  # geom_line(data = path_df %>% filter(Noise_sd %in% c(0.5, 0.8, 1.0)),
  #           mapping = aes(model_dim, rms, color = treatment, shape = ind_rep)) +
  geom_point(data = path_df %>% filter(Noise_sd %in% c(0.5, 0.8, 1.0)),
             mapping = aes(model_dim, rms, color = treatment, shape = as.factor(ind_rep))) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  ylab("Root mean square") + xlab("Model dimension") +
  scale_y_log10() +
  geom_point(data = path_df %>% group_by(treatment, Noise_sd, ind_rep) %>%
               slice(c(which.min(aic), which.min(bic), which.min(gcv))) %>%
               filter(Noise_sd %in% c(0.5, 0.8, 1.0)) %>% ungroup(),
             mapping = aes(model_dim, rms, colour = treatment),
             size = 4, fill = "white", shape = rep(23:25, 6))
ggsave(file = "simu/figure/dataset_pc_province/dataset_pc_province_rms_model_dim_crit.pdf",
       width = 6, height = 11)

ggplot(path_df %>% filter(Noise_sd %in% c(0.5, 0.8, 1.0)),
       aes(model_dim, rms, color = treatment)) +
  facet_grid(Noise_sd ~., labeller = labeller(Noise_sd = label_both)) +
  geom_line() + geom_point() +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = c(0.15, 0.10)) +
  ylab("Root mean square") + xlab("Model dimension") +
  scale_x_log10() + scale_y_log10()
ggsave(file = "simu/figure/dataset_pc_province/dataset_pc_province_rms_path_model_dim.pdf",
       width = 6, height = 11)

ind = 5
plot(dim_path_flsa[[ind]], rms_path_flsa[[ind]], type = "l", col = "green")
lines(dim_path[[ind]], rms_path[[ind]], type = "l", col = "red")

par(mfrow = c(1, 2))
plot(gamma_array[, ind, 1])
points(full_infer_res_flsa[[ind]]$result[which.min(full_infer_res[[ind]]$aic), ], col = "blue", pch = 3)
points(dataset$baseline, col = "red")
plot(gamma_array[, ind, 1])
points(full_infer_res[[ind]]$result[which.min(full_infer_res[[ind]]$aic), ], col = "green", pch = 3)
points(dataset$baseline, col = "red")

# run inference with multiple realizations -------------------------------------------
# infer_res <- apply(gamma_array, c(2, 3), infer, graph = graph, lambda = lambda,
#                      baseline = dataset$baseline)
cl <- makeCluster(getOption("cl.cores", detectCores() - 2))
parallel::clusterExport(cl = cl, varlist = ls())
parallel::clusterEvalQ(cl = cl, {
  library(tidyverse); library(gdata); library(RColorBrewer)
  library(rgeos); library(rgdal)
  library(glmnet); library(igraph); library(sf); library(Matrix)
  library(pryr); library(limSolve)
  library(scales); library(latex2exp); library(flsa)
  library(kableExtra)})
infer_res <- parApply(cl, gamma_array, c(2, 3), infer, graph = graph,
                      lambda = lambda, baseline = dataset$baseline)
infer_res_flsa <- parApply(cl, gamma_array, c(2, 3), infer_flsa, graph = graph,
                           lambda = lambda, baseline = dataset$baseline)
stopCluster(cl = cl)
# saveRDS(infer_res, file = paste0("simu_results/", region, "_", area, "_pc_", zone, "_infer_res.rds"))
# saveRDS(infer_res_flsa, file = paste0("simu_results/", region, "_", area, "_", zone, "_infer_res_flsa.rds"))
infer_res <- readRDS(file = paste0("simu/simu_results/", region, "_", area, "_", zone, "_infer_res.rds"))
infer_res_flsa <- readRDS(file = paste0("simu/simu_results/", region, "_", area, "_", zone, "_infer_res_flsa.rds"))


# Post-treatment
## Root mean square
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

## Dimension
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
dim_rms_df <- inner_join(dim_df, rms_df, by = c("sigma_sq", "criterion", "ind_rep", "treatment"))

# save simulation results -------------------------------------------
# saveRDS(dim_rms_df, paste0("simu_results/infer_", region, "_", area, "_pc_", zone, "dim_rms.rds"))
dim_rms_df <- readRDS(file = paste0("simu/simu_results/infer_", region, "_", area, "_pc_", zone, "dim_rms.rds"))

### dimension table -------------------------------------------
rms0 <- sapply(seq_along(sigma_sq), function(a)
  sqrt(sum((gamma_array[, a, 1] - dataset$baseline) ^ 2 ) / length(gamma))
)

ccc <- dim_rms_df %>% 
  filter(criterion == "bic") %>%
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
saveRDS(eee, file = "simu/table/dataset_pc_province/dataset_pc_province_rms_dim_df.rds")
saveRDS(ddd, file = "simu/table/dataset_pc_province/dataset_pc_province_rms_dim_kable.rds")
ddd %>% cat(file = "simu/table/dataset_pc_province/dataset_pc_province_rms_dim_table.tex")

### rms table -------------------------------------------
rms0 <- sapply(seq_along(sigma_sq), function(a)
  sqrt(sum((gamma_array[, a, 1] - dataset$baseline) ^ 2 ) / length(gamma))
)
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

bbb %>% cat(file = "simu/table/dataset/pc/dataset_pc_rms_table.tex")

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


## Run agraph once and plot results

for (ind_sigma in seq_along(sigma_sq)) {
  # Perform inference once
  res <- agraph(gamma = gamma_array[, ind_sigma, 1],
                graph = graph,
                lambda = lambda,
                weights = NULL,
                shrinkage = TRUE,
                tol = 1e-8,
                itermax = 50000)
  connlist <- graph2connlist(graph)
  res2 <- flsa(gamma_array[, ind_sigma, 1], lambda2 = lambda, connListObj = connlist)
  res3_obj <- flsa(gamma_array[, ind_sigma, 1], lambda2 = NULL, connListObj = connlist)
  res3 <- flsaGetSolution(res3_obj, lambda2 = lambda)
  matplot(res3, type = "l")
  plot(apply(res3, 1, function(a) length(unique(a))))
  
  dataset$theta <- res$result[which.min(res$bic), ]
  dataset$theta_flsa <- res2[which.min(res$bic), ]
  
  
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
  ggsave(file = paste0("simu/figure/dataset/pc/utr_neigh_pc_reg_path_sigma",
                       ind_sigma,".pdf"), width = width_pdf, height = height_pdf)
  
  ggplot(data = res_path_flsa, aes(lambda, value, group = area)) +
    geom_line(size = 0.2) +
    scale_x_log10() +
    xlab(TeX("$\\log_{10}(\\lambda)$")) + ylab("Parameters") +
    theme_minimal() +
    geom_vline(xintercept = lambda[sapply(res[c("aic", "bic", "gcv")], which.min)], col = "red") +
    theme(legend.position = "none")
  ggsave(file = paste0("simu/figure/dataset/pc/utr_neigh_pc_reg_path_flsa_sigma",
                       ind_sigma,".pdf"), width = 6.5, height = height_pdf)
  
  # matplot(log10(lambda), res$result, type = "l", lwd = 0.3, lty = 1,
  #         col = "black", ylab = "Parameter values", xlab = TeX("$\\log_{10}(\\lambda)$"))
  # abline(v = log10(lambda)[sapply(res[c("aic", "bic", "gcv")], which.min)],
  #        col = "red")
  
  dataset_fused <- dataset %>%
    mutate(group = cut(theta, breaks = sort(unique(theta)),
                       include.lowest = TRUE)) %>%
    group_by(group) %>%
    dplyr::summarize(theta = mean(theta), do_union = TRUE)
  dataset_baseline_fused <- dataset %>%
    mutate(group = cut(baseline, breaks = sort(unique(baseline)),
                       include.lowest = TRUE)) %>%
    group_by(group) %>%
    dplyr::summarize(baseline = mean(baseline), do_union = TRUE)
  
  # Plot segmentations without color
  # op <- par(mar = rep(0, 4), mfrow = c(1, 2))
  # plot(st_geometry(dataset_fused), main = "", key.pos = NULL)
  # plot(st_geometry(dataset), add = TRUE, lwd = 0.2, border = "grey")
  # plot(st_geometry(dataset_baseline_fused), main = "", key.pos = NULL)
  # plot(st_geometry(dataset), add = TRUE, lwd = 0.2, border = "grey")
  # par(op)
  # dev.copy2pdf(file = paste0("figure/dataset/pc/utr_neigh_pc_baseline_vs_seg_sigma",
  #                            ind_sigma,".pdf"), width = 7, height = 5)
  
  # Plot segmentation overlay
  op <- par(mar = rep(0, 4), mfrow = c(1, 1))
  plot(st_geometry(dataset_fused), main = "", key.pos = NULL,
       border = alpha("red", 0.8))
  plot(st_geometry(dataset_baseline_fused), main = "", key.pos = NULL,
       border = alpha("black", 0.8), add = TRUE)
  plot(st_geometry(dataset), add = TRUE, lwd = 0.2, border = "grey")
  par(op)
  dev.copy2pdf(file = paste0("simu/figure/dataset/pc/utr_neigh_pc_baseline_seg_overlay_sigma",
                             ind_sigma,".pdf"), width = width_pdf, height = height_pdf)
  
  # Define palette of initial signal and baseline signal
  x.breaks <- seq(from = min(gamma_array[, ind_sigma, 1]),
                  to = max(gamma_array[, ind_sigma, 1]),
                  by = 0.01)
  palette_signal <- pal_div(n = length(x.breaks) - 1)[cut(x = gamma_array[, ind_sigma, 1],
                                                          breaks = x.breaks, right = FALSE,
                                                          include.lowest = TRUE, labels = FALSE)]
  palette_baseline_fused <- pal_div(n = length(x.breaks) - 1)[cut(x = dataset_baseline_fused$baseline,
                                                                  breaks = x.breaks, right = FALSE,
                                                                  include.lowest = TRUE, labels = FALSE)]
  # Define palette of segmented signal
  palette_seg <- pal_div(n = length(x.breaks) - 1)[cut(x = dataset$theta,
                                                       breaks = x.breaks, right = FALSE,
                                                       include.lowest = TRUE, labels = FALSE)]
  palette_seg_fused <- pal_div(n = length(x.breaks) - 1)[cut(x = dataset_fused$theta,
                                                             breaks = x.breaks, right = FALSE,
                                                             include.lowest = TRUE, labels = FALSE)]
  
  # Plot unsegmented signal
  dataset$gamma <- gamma_array[, ind_sigma, 1]
  op <- par(mar = rep(0, 4))
  plot(dataset["gamma"], col = palette_signal,
       main = "", lwd = 0.2, border = alpha("black", 0.6))
  par(op)
  dev.copy2pdf(file = paste0("simu/figure/dataset/pc/utr_neigh_pc_unseg_sigma",
                             ind_sigma,".pdf"), width = width_pdf, height = height_pdf)
  
  # Plot baseline signal
  op <- par(mar = rep(0, 4))
  plot(st_geometry(dataset_baseline_fused), main = "",
       lwd = 1.2, border = "black", col = palette_baseline_fused)
  plot(st_geometry(dataset), col = NA, add = TRUE,
       main = "", lwd = 0.2, border = alpha("black", 0.6))
  par(op)
  dev.copy2pdf(file = paste0("simu/figure/dataset/pc/utr_neigh_pc_baseline_sigma",
                             ind_sigma,".pdf"), width = width_pdf, height = height_pdf)
  
  # Plot segmented signal
  op <- par(mar = rep(0, 4))
  plot(st_geometry(dataset_fused), lwd = 1.2, border = "black",
       col = palette_seg_fused)
  plot(st_geometry(dataset), add = TRUE, lwd = 0.2, border = alpha("black", 0.6))
  par(op)
  dev.copy2pdf(file = paste0("simu/figure/dataset/pc/utr_neigh_pc_seg_sigma",
                             ind_sigma,".pdf"), width = width_pdf, height = height_pdf)
}


# Visual inspection: what is the signal-to-noise ratio?
plot(dataset$baseline)
points(gamma_array[,5,1], col = "red", pch = 3)
plot(dataset$baseline)
points(gamma_array[,8,1], col = "green", pch = 3)

