library(tidyverse); library(gdata); library(RColorBrewer)
library(glmnet); library(igraph); library(sf); library(Matrix)
library(pryr); library(limSolve); library(graphseg)
library(scales); library(latex2exp); library(flsa)
library(kableExtra)

# Load data ---------------------------------------------------------------
tidy_df <- bind_rows(
  readRDS(file = "simu/simu_results/neth_municip_pc_province_tidy_df.rds") %>% 
    mutate(dataset = "dataset_1"),
  readRDS(file = "simu/simu_results/utrecht_district_pc_municip_tidy_df.rds") %>% 
    mutate(dataset = "dataset_2"),
  readRDS(file = "simu/simu_results/utrecht_neigh_pc_province_tidy_df.rds") %>% 
    mutate(dataset = "dataset_3"),
  readRDS(file = "simu/simu_results/utrecht_neigh_pc_municip_tidy_df.rds") %>% 
    mutate(dataset = "dataset_4"),
  readRDS(file = "simu/simu_results/utrecht_neigh_pc_district_tidy_df.rds") %>% 
    mutate(dataset = "dataset_5"),
  readRDS(file = "simu/simu_results/neth_neigh_pc_municip_tidy_df.rds") %>% 
    mutate(dataset = "dataset_6")
)

rms0_df <- tidy_df %>% select(sigma_sq, dataset, rms0) %>% 
  distinct()

# Table with noise sd = 0.7 -----------------------------------------------
# temp <- tidy_df %>% filter(sigma_sq == 0.7) %>% 
#   select(-rms0) %>% 
#   filter(criterion == "aic") %>% filter(dataset == "dataset_3") %>% 
#   group_by(treatment) %>% 
#   summarize(m_rms = mean(rms), m_clust_score = mean(clust_score), m_dim = mean(dim))
# temp
# 
# temp2 <- tidy_df %>% filter(sigma_sq == 0.1, dataset == "dataset_6", 
#                             criterion == "aic") %>% 
#   group_by(treatment) %>% summarise(mean_dim = mean(dim))
# temp2

a <- tidy_df %>% filter(sigma_sq == 0.7) %>% 
  filter(criterion == "aic") %>% 
  select(-rms0, -sigma_sq, -criterion)
b <- a %>% gather(key = "quantity", value = "value", -c(ind_rep, treatment, dataset)) %>%
  group_by(treatment, dataset, quantity)
c <- b %>%
  dplyr::summarize(mean = mean(value),
                   sd = sd(value), .groups = "keep") %>%
  mutate(mean = signif(mean, digits = 3) %>% sprintf("%.3f", .),
         sd = signif(sd, digits = 2) %>% sprintf("%.2f", .)) %>%
  mutate(value_err = paste('$', mean, "\\pm", sd, '$')) %>% 
  ungroup()

d <- c %>% 
  mutate(column = paste0(quantity, "_", treatment)) %>% 
  select(-mean, -sd, -treatment, -quantity) %>% 
  spread(key = "column", value = "value_err")

min_col <- rbind(c %>% select(-c(value_err, sd)) %>%
                           reshape2::acast(dataset ~ quantity + treatment) %>%
                           apply(1, function(a) which.max(a[1:2])),
                         c %>% select(-c(value_err, sd)) %>%
                           reshape2::acast(dataset ~ quantity + treatment) %>%
                           apply(1, function(a) which.min(a[3:4]) + 2),
                         c %>% select(-c(value_err, sd)) %>%
                           reshape2::acast(dataset ~ quantity + treatment) %>%
                           apply(1, function(a) which.min(a[5:6]) + 4))
for (ind1 in 1:ncol(min_col)) {
  for (ind2 in 1:nrow(min_col)) {
    d[ind1, min_col[ind2, ind1] + 1] <-
      paste0("\\boldmath{",
             d[ind1, min_col[ind2, ind1] + 1],
             "}")
  }
}
# prod_table_rms_dim_clust_score <- c %>% 
#   reshape2::acast(dataset ~ quantity + treatment) %>%
#   as.data.frame(stringsAsFactors = FALSE)

e <- d %>% mutate(p = c(390, 650, 2955, 2955, 2955, 12920),
                  true_dim = c(12, 99, 6, 99, 650, 3068)) %>% 
  left_join(rms0_df %>% filter(sigma_sq == 0.7), by = "dataset") %>% 
  mutate(rms0 = signif(rms0, digits = 2)) %>% 
  select(-dataset, -sigma_sq) %>% 
  relocate(p, dim_aridge, dim_flsa, true_dim, rms_aridge, rms_flsa, rms0,
           clust_score_aridge, clust_score_flsa) %>% 
  mutate(p = paste0("Dataset ", 1:6, "\n(", p, ")")) %>% 
  mutate(p = kableExtra::linebreak(p, align = "r"))

rms_dim_clust_score_kable <- e %>% kable(
  format = "latex", label = "rms_dim_clust_score_table", caption = 
    "Model dimension (lower is better, in bold), RMSE (lower is better, in bold), and adjusted Rand index (higher is better, in bold) of the adaptive ridge compared to Hoefling's fused lasso signal approximator (lfsa), using the AIC. The reported values are the means over 100 repetitions, $\\pm$ the standard error. The noise standard deviation is $\\sigma = 0.7$.",
  booktabs = TRUE, escape = FALSE, linesep = "",
  col.names = linebreak(c("Dataset \n (p)",
                          c("Adaptive \n Ridge", "FLSA"), "True \n model \n dim (q)",
                          c("Adaptive \n Ridge", "FLSA"), "No \n Regula- \n rization",
                          c("Adaptive \n Ridge", "FLSA")), align = "c")) %>% 
  add_header_above(c(" ","Model \n dimension" = 3, "RMSE" = 3, "Adjusted Rand \n index" = 2)) %>% 
  kable_styling(latex_options = c("striped", "repeat_header", "scale_down"))

cat(rms_dim_clust_score_kable, file = "../paper/submission_csda/revised_manuscript/rms_dim_clust_score_table.tex")
cat(rms_dim_clust_score_kable, file = "../paper/jan_vivien_segmentation/rms_dim_clust_score_table.tex")
cat(rms_dim_clust_score_kable, file = "simu/table/rms_dim_clust_score_table.tex")


# Table with all noise sd ------------------------------------------

a2 <- tidy_df %>% 
  filter(criterion == "aic") %>% 
  filter(sigma_sq %in% c(0.1, 0.3, 0.7, 0.9, 1.1, 2, 5)) %>% 
  select(-rms0, -criterion)
b2 <- a2 %>% gather(key = "quantity", value = "value", -c(sigma_sq, ind_rep, treatment, dataset)) %>%
  group_by(treatment, sigma_sq, dataset, quantity)
c2 <- b2 %>%
  dplyr::summarize(mean = mean(value),
                   sd = sd(value), .groups = "keep") %>%
  mutate(mean = signif(mean, digits = 3) %>% as.character,
         sd = signif(sd, digits = 2) %>% as.character) %>%
  mutate(value_err = paste('$', mean, "\\pm", sd, '$')) %>% 
  ungroup()
d2 <- c2 %>% 
  mutate(column = paste0(quantity, "_", treatment)) %>% 
  select(-mean, -sd, -treatment, -quantity) %>% 
  spread(key = "column", value = "value_err")

min_col2 <- rbind(c2 %>% select(-c(value_err, sd)) %>%
                   reshape2::acast(sigma_sq + dataset ~ quantity + treatment) %>%
                   apply(1, function(a) which.max(a[1:2])), # clust_score
                 c2 %>% select(-c(value_err, sd)) %>%
                   reshape2::acast(sigma_sq + dataset ~ quantity + treatment) %>%
                   apply(1, function(a) which.min(a[3:4]) + 2), # model_dim
                 c2 %>% select(-c(value_err, sd)) %>%
                   reshape2::acast(sigma_sq + dataset ~ quantity + treatment) %>%
                   apply(1, function(a) which.min(a[5:6]) + 4)) # rmse

for (ind1 in 1:ncol(min_col2)) {
  for (ind2 in 1:nrow(min_col2)) {
    d2[ind1, min_col2[ind2, ind1] + 2] <-
      paste0("\\boldmath{",
             d2[ind1, min_col2[ind2, ind1] + 2],
             "}")
  }
}
d2

p_q_df <- tibble(dataset = paste0("dataset_", 1:6),
                 p = c(390, 650, 2955, 2955, 2955, 12920),
                 true_dim = c(12, 99, 6, 99, 650, 3068))

e2 <- d2 %>% left_join(p_q_df, by = "dataset") %>% 
  left_join(rms0_df) %>% 
  mutate(rms0 = signif(rms0, digits = 2)) %>% 
  filter(sigma_sq %in% c(0.1, 0.3, 0.7, 0.9, 1.1, 2, 5)) %>% 
  select(-dataset, -sigma_sq) %>% 
  relocate(p, dim_aridge, dim_flsa, true_dim, rms_aridge, rms_flsa, rms0,
           clust_score_aridge, clust_score_flsa) %>% 
  mutate(p = paste0("Dataset ", 1:6, "\n(", p, ")")) %>% 
  mutate(p = kableExtra::linebreak(p, align = "r"))

rms_dim_clust_score_kable2 <- e2 %>% kable(
  format = "latex", label = "rms_dim_clust_score_table_whole", caption = 
    "Same table as Table 2 for varying values of $\\sigma$. Model dimension, RMSE, and adjusted Rand index (higher is better) of the adaptive ridge and the lfsa, using the AIC. The reported values are the means over 100 repetitions, $\\pm$ the standard error.",
  booktabs = TRUE, escape = FALSE, linesep = "",
  longtable = TRUE,
  col.names = linebreak(c("Dataset \n (p)",
                          c("Adaptive \n Ridge", "FLSA"), "True \n model \n dim (q)",
                          c("Adaptive \n Ridge", "FLSA"), "No \n Regula- \n rization",
                          c("Adaptive \n Ridge", "FLSA")), align = "c")) %>% 
  add_header_above(c(" ","Model \n dimension" = 3, "RMSE" = 3, "Adjusted Rand \n index" = 2)) %>% 
  kable_styling(latex_options = c("striped", "repeat_header"), font_size = 6) %>% 
  pack_rows("$\\sigma = 0.1$", 1, 6, hline_after = TRUE, hline_before = F, escape = FALSE,latex_gap_space = "0.3em") %>%
  pack_rows("$\\sigma = 0.3$", 7, 12, hline_after = TRUE, hline_before = F, escape = FALSE,latex_gap_space = "0.3em") %>%
  pack_rows("$\\sigma = 0.7$", 13, 18, hline_after = TRUE, hline_before = F, escape = FALSE, latex_gap_space = "0.3em") %>%
  pack_rows("$\\sigma = 0.9$", 19, 24, hline_after = TRUE, hline_before = F, escape = FALSE, latex_gap_space = "0.3em") %>%
  pack_rows("$\\sigma = 1.1$", 25, 30, hline_after = TRUE, hline_before = F, escape = FALSE, latex_gap_space = "0.3em") %>%
  pack_rows("$\\sigma = 2$", 31, 36, hline_after = TRUE, hline_before = F, escape = FALSE, latex_gap_space = "0.3em") %>%
  pack_rows("$\\sigma = 5$", 37, 42, hline_after = TRUE, hline_before = F, escape = FALSE, latex_gap_space = "0.3em")
  
cat(rms_dim_clust_score_kable2, file = "../paper/submission_csda/revised_manuscript/rms_dim_clust_score_table_whole.tex")
system("ex ../paper/submission_csda/revised_manuscript/rms_dim_clust_score_table_whole.tex <<eof
2 insert \n \\setlength\\LTleft{-1.6cm}\n.\nxit\neof")
cat(rms_dim_clust_score_kable2, file = "../paper/jan_vivien_segmentation/rms_dim_clust_score_table_whole.tex")
system("ex ../paper/jan_vivien_segmentation/rms_dim_clust_score_table_whole.tex <<eof
2 insert \n \\setlength\\LTleft{-1.6cm}\n.\nxit\neof")
cat(rms_dim_clust_score_kable2, file = "simu/table/rms_dim_clust_score_table_whole.tex")
system("ex simu/table/rms_dim_clust_score_table_whole.tex <<eof
2 insert \n \\setlength\\LTleft{-1.6cm}\n.\nxit\neof")
