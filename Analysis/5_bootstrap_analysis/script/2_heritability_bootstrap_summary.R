library(lme4)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# source the heritability calculation function
source("../../../Code/heritability_fun.R")
source("../../../Code/plot_fun.R")

outdir = "../result/"

#-------read h2 results from 6_heritability_bootstrap.R----------
h2_results_long <- read.csv(paste0(outdir, "h2_bootstrap.csv"))
results_summary <- h2_results_long %>%
  group_by(n_geno, n_rep, metric) %>%
  summarise(
    h2_mean = mean(h2, na.rm = TRUE),
    h2_sd   = sd(h2,   na.rm = TRUE),
    h2_low_q  = quantile(h2, 0.025, na.rm = TRUE),
    h2_high_q = quantile(h2, 0.975, na.rm = TRUE)
  )



summary_list <- results_summary %>% 
  group_by(metric) %>% 
  group_split()

main_metrics_order <- unique(results_summary$metric)
names(summary_list) <- main_metrics_order# not sure if in the right order



# ------Plot H² vs #genotypes for each metric------

for (met in names(summary_list)) {
  df_met <- summary_list[[met]]
  
  p <- ggplot(df_met,
              aes(x = n_geno,
                  y = h2_mean,
                  color = factor(n_rep),
                  group = n_rep)) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = h2_low_q, ymax = h2_high_q),
                  width = 1, alpha = 0.6) +
    labs(
      title = paste("Heritability vs sample size for", met),
      x = "Number of genotypes",
      y = "Heritability (H²)",
      color = "Replicates / genotype"
    ) +
    theme_bw()
  
  print(p)  # show on screen
  
  # 👇 optional: save to file
  ggsave(
    filename = paste0(outdir,"/H2_bootstrap_summary_", met, ".png"),
    plot = p,
    width = 6, height = 4, dpi = 300
  )
}





# ------Plot H²'s surface plot  for all metrics------
# 1. Start from the full data
# h2_results_long has: metric, n_geno, n_rep, h2

#### use function to plot
plot_surface_z(results_summary,z_value='h2_sd',mean_col='h2_mean')


df_summary <- h2_results_long %>%
  group_by(metric, n_geno, n_rep) %>%
  summarise(
    mean_h2 = mean(h2, na.rm = TRUE),
    'Negtive Variance'  = 1-var(h2,  na.rm = TRUE),
    IQR  = IQR(h2, na.rm = TRUE),# it's not sd, but IQR
    Spread  = quantile(h2, 0.975) - quantile(h2, 0.025),
    .groups = "drop"
  )

plot_surface_z(df_summary,z_value='Negtive Variance',mean_col='mean_h2',
               outdir=outdir,
               file_name = "H2_surface_plot_var.html"
               )

plot_surface_z(df_summary,z_value='Spread',mean_col='mean_h2')



