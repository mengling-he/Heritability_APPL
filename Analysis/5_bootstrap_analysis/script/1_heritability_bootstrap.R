library(lme4)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# source the heritability calculation function
source("../../../Code/heritability_fun.R")

outdir = "../result/"


# ----read data----
##### read data ####
features_df0_top_side <- read.csv("../data/side_top_features_vit_linearimputed_day.csv",stringsAsFactors = FALSE)
colnames(features_df0_top_side)

##### delete dead  #####
slow_growers <- strsplit(trimws(readLines("../data/slow_grower.txt")), ",")[[1]]
features_df_top_side<- features_df0_top_side  %>% filter(! Plant.Info %in% slow_growers)

##### generate a Day column  #####
features_df_top_side$Date <- as.Date(features_df_top_side$Date)
features_df_top_side$Day <- as.numeric(features_df_top_side$Date - min(features_df_top_side$Date))
colnames(features_df_top_side)

main_metrics <- c("side_height","side_width","side_solidity","side_area","top_blue_yellow_mean","top_saturation_mean")
columns_to_kept <- c("Plant.Info","Date","Genotype","Replicate","Day",main_metrics)
df_final <- features_df_top_side[,columns_to_kept]


# ----design settings----
genotypes_list  <- c(10, 20,30,40,50,60)
replicates_list <- c(2,3,4,5, 6)
max_reps        <- max(replicates_list)   # = 6
n_repeat        <- 100               # number of random rounds

##### replicate order for genotypes #####
all_genotypes <- unique(df_final$Genotype)


##### Pre-compute a replicate order per genotype (for nested reps)#####
rep_pool_df <- df_final %>%
  group_by(Genotype) %>%
  summarise(
    rep_ids = list(
      sample(
        unique(Replicate),
        size = min(max_reps, dplyr::n_distinct(Replicate))
      )
    ),
   .groups = "drop"
  )

rep_pool_list <- setNames(rep_pool_df$rep_ids,
                          as.character(rep_pool_df$Genotype))


######## calculate H2 for nested genotypes + nested replicates: repeat 100 times  #######

set.seed(999)   # global seed for reproducibility

# for 100 repeats, 2 seed streams
seed_geno_list <- sample(1e8, n_repeat)
seed_reps_list <- sample(1e8, n_repeat)

res_list <- vector("list", length = n_repeat)

for (r in seq_len(n_repeat)) {
  ## ---- 1. Random order of all genotypes for this repeat ----
  set.seed(seed_geno_list[r])
  geno_perm <- sample(all_genotypes, length(all_genotypes))
  
  ## ---- 2. Use a separate seed stream for replicate sampling in this repeat ----
  set.seed(seed_reps_list[r])
  
  round_results <- list()
  
  for (n_geno in genotypes_list) {
    geno_subset <- geno_perm[1:n_geno]
    
    for (n_rep in replicates_list) {
      
      # For this (n_geno, n_rep), sample replicates *randomly* for each genotype
      df_sub <- df_final %>%
        dplyr::filter(Genotype %in% geno_subset) %>%
        dplyr::group_by(Genotype) %>%
        dplyr::group_modify(~ {
          # .x is the data for this genotype
          available_reps <- unique(.x$Replicate)
          
          # how many to take for this genotype
          n_take <- min(length(available_reps), n_rep)
          
          # RANDOM choice of replicates, independent for each n_rep
          reps_keep <- sample(available_reps, n_take, replace = FALSE)
          
          dplyr::filter(.x, Replicate %in% reps_keep)
        }) %>%
        dplyr::ungroup()
      
      if (length(unique(df_sub$Genotype)) < 2) {
        h2_df_round <- data.frame(
          round  = r,
          n_geno = n_geno,
          n_rep  = n_rep,
          metric = main_metrics,
          h2     = NA_real_
        )
      } else {
        h2_results_df <- calc_H2_fixedtime(df_sub, main_metrics)
        
        h2_df_round <- h2_results_df %>%
          dplyr::mutate(
            round  = r,
            n_geno = n_geno,
            n_rep  = n_rep
          ) %>%
          dplyr::rename(
            metric = Trait,
            h2     = H2
          )
      }
      
      round_results[[length(round_results) + 1]] <- h2_df_round
    }
  }
  
  res_list[[r]] <- dplyr::bind_rows(round_results)
}


# ------finalize-------
results_long <- bind_rows(res_list)

write.csv(
  results_long,
  file = paste0(outdir, "h2_bootstrap.csv"),
  row.names = FALSE
)






# ----h2 results summary----
# refer to 6_heritability_bootstrap_summary.R



# --- some other exploration ----
target_metric <- "top_blue_yellow_mean"  
#for (target_metric in main_metrics){
df_met <- h2_results_long %>%
  filter(metric == target_metric, !is.na(h2))

# 
# df_met_lm <- df_met %>%
#   mutate(
#     n_geno_f = factor(n_geno),
#     n_rep_f  = factor(n_rep)
#   )
# 
# fit <- lm(h2 ~ n_geno_f * n_rep_f, data = df_met_lm)
# 
# anova_tab <- anova(fit)
# anova_tab
# 
# anova_var <- anova_tab %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column("source") %>%
#   mutate(
#     prop_SS = `Sum Sq` / sum(`Sum Sq`, na.rm = TRUE)
#   )
# 
# anova_var
# ggplot(anova_var %>% filter(source != "Residuals"),
#        aes(x = source, y = prop_SS)) +
#   geom_col() +
#   coord_flip() +
#   labs(
#     title = paste("Variance decomposition of H² for", target_metric),
#     x = "Source",
#     y = "Proportion of variance in H² explained"
#   ) +
#   theme_bw()
# 

summary_met <- df_met %>%
  group_by(n_geno, n_rep) %>%
  summarise(
    h2_mean  = mean(h2, na.rm = TRUE),
    h2_sd    = sd(h2,   na.rm = TRUE),
    h2_low_q = quantile(h2, 0.025, na.rm = TRUE),
    h2_high_q= quantile(h2, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(summary_met,
       aes(x = n_geno,
           y = h2_mean,
           color = factor(n_rep),
           group = n_rep)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = h2_low_q, ymax = h2_high_q),
                width = 1, alpha = 0.6) +
  labs(
    title = paste("Sensitivity of H² to design for", target_metric),
    x = "Number of genotypes",
    y = "Heritability (H²)",
    color = "Replicates / genotype"
  ) +
  theme_bw()


ggplot(summary_met,
       aes(x = n_geno,
           y = n_rep,
           fill = h2_mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(h2_mean, 2)), size = 3) +
  scale_y_continuous(breaks = sort(unique(summary_met$n_rep))) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = paste("H² design surface for", target_metric),
    x = "Number of genotypes",
    y = "Replicates / genotype",
    fill = "Mean H²"
  ) +
  theme_bw()


summary_met2 <- summary_met %>%
  mutate(ci_width = h2_high_q - h2_low_q)

ggplot(summary_met2,
       aes(x = n_geno,
           y = n_rep,
           fill = ci_width)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(ci_width, 2)), size = 3) +
  scale_y_continuous(breaks = sort(unique(summary_met2$n_rep))) +
  scale_fill_gradient(low = "white", high = "darkred") +
  labs(
    title = paste("H² CI width surface for", target_metric),
    x = "Number of genotypes",
    y = "Replicates / genotype",
    fill = "CI width"
  ) +
  theme_bw()

