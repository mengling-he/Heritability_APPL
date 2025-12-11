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


df_summary <- h2_results_long %>%
  group_by(metric, n_geno, n_rep) %>%
  summarise(
    mean_h2 = mean(h2, na.rm = TRUE),
    neg_var_h2  = 1-var(h2,  na.rm = TRUE),
    #neg_sd_h2  = IQR(h2, na.rm = TRUE),# it's not sd, but IQR
    neg_sd_h2  = quantile(h2, 0.975) - quantile(h2, 0.025),
    .groups = "drop"
  )

# Axes (assume same grid for all metrics)
x_vals <- sort(unique(df_summary$n_geno))  # n_geno
y_vals <- sort(unique(df_summary$n_rep))   # n_rep

metrics <- sort(unique(df_summary$metric))

# 2. Precompute matrices for each metric
var_mats  <- list()
sd_mats  <- list()
mean_mats <- list()
highlight_list <- list()

for (m in metrics) {
  m_df <- df_summary %>% filter(metric == m)
  
  # matrices: rows = y_vals (n_rep), cols = x_vals (n_geno)
  var_mat  <- matrix(NA_real_, nrow = length(y_vals), ncol = length(x_vals))
  sd_mat  <- matrix(NA_real_, nrow = length(y_vals), ncol = length(x_vals))
  mean_mat <- matrix(NA_real_, nrow = length(y_vals), ncol = length(x_vals))
  
  for (i in seq_along(x_vals)) {
    for (j in seq_along(y_vals)) {
      sub_df <- m_df %>%
        filter(n_geno == x_vals[i], n_rep == y_vals[j])
      if (nrow(sub_df) == 1) {
        var_mat[j, i]  <- sub_df$neg_var_h2
        sd_mat[j, i]  <- sub_df$neg_sd_h2
        mean_mat[j, i] <- sub_df$mean_h2
      }
    }
  }
  
  var_mats[[m]]  <- var_mat
  sd_mats[[m]]  <- sd_mat
  mean_mats[[m]] <- mean_mat
  
  # highlight path on this surface
  highlight_list[[m]] <- data.frame(
    n_geno = c(20, 30, 40, 60),
    n_rep  = c(6,  4,  3,  2)
  ) %>%
    left_join(m_df, by = c("n_geno", "n_rep"))
}

# 3. Build the plot with a dropdown
p <- plot_ly()

for (k in seq_along(metrics)) {
  m <- metrics[k]
  var_mat  <- var_mats[[m]]
  sd_mat  <- sd_mats[[m]]
  mean_mat <- mean_mats[[m]]
  hl       <- highlight_list[[m]]
  
  visible_init <- (k == 1)  # only first metric visible at start
  
  # Surface: z = variance, color = mean
  p <- p %>%
    add_surface(
      x = x_vals,
      y = y_vals,
      z = sd_mat,
      surfacecolor = mean_mat,
      colorscale = "Viridis",
      showscale = (k == 1),   # show colorbar only for active metric
      colorbar = list(title = "Mean H²"),
      visible = visible_init,
      name = paste(m, "surface")
    )
  
  # Line + markers on the surface
  p <- p %>%
    add_trace(
      data = hl,
      x = ~n_geno,
      y = ~n_rep,
      z = ~neg_sd_h2,
      type = "scatter3d",
      mode = "lines+markers",
      marker = list(size = 6, color = "red"),
      line   = list(color = "red", width = 5),
      visible = visible_init,
      showlegend = FALSE,
      inherit = FALSE
    )
  
  # Text labels: mean(H²)
  p <- p %>%
    add_text(
      data = hl,
      x = ~n_geno,
      y = ~n_rep,
      z = ~neg_sd_h2,
      text = ~sprintf("mean H2=%.3f", mean_h2),
      type = "scatter3d",
      mode = "text",
      visible = visible_init,
      showlegend = FALSE,
      inherit = FALSE
    )
}

# Build dropdown buttons
n_metrics <- length(metrics)
n_traces  <- 3 * n_metrics  # surface + line + text per metric

buttons <- lapply(seq_len(n_metrics), function(k) {
  vis <- rep(FALSE, n_traces)
  idx <- ((k - 1) * 3 + 1):(k * 3)  # traces for metric k
  vis[idx] <- TRUE
  
  list(
    method = "update",
    args = list(
      list(visible = vis),
      list(title = paste("Spread Surface of H² -", metrics[k]))
    ),
    label = metrics[k]
  )
})

p <- p %>%
  layout(
    title = paste("Spread Surface of H² -", metrics[1]),
    scene = list(
      xaxis = list(title = "n_geno"),
      yaxis = list(title = "n_rep"),
      #zaxis = list(title = "1-Var(H²)")
      zaxis = list(title = "CI spread(H²)")
    ),
    updatemenus = list(
      list(
        type = "dropdown",
        x = 0.05,
        y = 1.05,
        showactive = TRUE,
        buttons = buttons
      )
    )
  )

p

htmlwidgets::saveWidget(p, paste0(outdir, "/H2_surface_plot_spread.html"))






# ------Plot H² vs #genotypes of heatmap for one metri: side_solidity------
df_trait <- h2_results_long %>% filter(metric == "side_solidity")

df_trait_sub <- df_trait[df_trait$n_geno==10&df_trait$n_rep==6,]


mean_h2_solidity <- mean(df_trait_sub$h2)
df_trait_sub$bias <- (df_trait_sub$h2-mean_h2_solidity)^2
1-sum(df_trait_sub$bias)/100
##### surface plot of variance #######
df_summary_trait <- df_trait %>%
  group_by(n_geno, n_rep) %>%
  summarise(
    mean_h2 = mean(h2, na.rm = TRUE),
    var_h2  = var(h2,  na.rm = TRUE),
    .groups = "drop"
  )

# unique sorted axes
x_vals <- sort(unique(df_summary_trait$n_geno))  # n_geno
y_vals <- sort(unique(df_summary_trait$n_rep))   # n_rep

var_mat  <- matrix(NA_real_, nrow = length(y_vals), ncol = length(x_vals))
mean_mat <- matrix(NA_real_, nrow = length(y_vals), ncol = length(x_vals))

for (i in seq_along(x_vals)) {
  for (j in seq_along(y_vals)) {
    sub_df <- df_summary_trait %>%
      filter(n_geno == x_vals[i], n_rep == y_vals[j])
    if (nrow(sub_df) == 1) {
      var_mat[j, i]  <- sub_df$var_h2
      mean_mat[j, i] <- sub_df$mean_h2
    }
  }
}

# # plot without path
# plot_ly(
#   x = x_vals,          # n_geno
#   y = y_vals,          # n_rep
#   z = var_mat,         # variance on z-axis
#   type = "surface",
#   surfacecolor = mean_mat,  # color encodes mean H²
#   colorscale = "Viridis",
#   colorbar = list(title = "Mean H²")
# ) %>%
#   layout(
#     title = "Variance Surface of H² (Color = Mean H²)",
#     scene = list(
#       xaxis = list(title = "n_geno"),
#       yaxis = list(title = "n_rep"),
#       zaxis = list(title = "Var(H²)")
#     )
#   )


# specific points to highlight
highlight_pts <- data.frame(
  n_geno = c(20, 30, 40, 60),
  n_rep  = c(6,  4,  3,  2)
) %>%
  left_join(df_summary, by = c("n_geno", "n_rep"))


p <- plot_ly(
  x = x_vals,
  y = y_vals,
  z = var_mat,
  type = "surface",
  surfacecolor = mean_mat,         # color encodes mean(H²)
  colorscale = "Viridis",
  colorbar = list(title = "Mean H²")
) %>%
  # line + markers on the surface
  add_trace(
    data = highlight_pts,
    x = ~n_geno,
    y = ~n_rep,
    z = ~var_h2,
    type = "scatter3d",
    mode = "lines+markers",
    marker = list(size = 6, color = "red"),
    line   = list(color = "red", width = 5),
    name = "Selected path",
    inherit = FALSE   # don't inherit surface attributes like colorscale/colorbar
  ) %>%
  # text labels showing mean(H²)
  add_text(
    data = highlight_pts,
    x = ~n_geno,
    y = ~n_rep,
    z = ~var_h2,
    text = ~sprintf("mean=%.3f", mean_h2),
    textposition = "top center",
    showlegend = FALSE,
    inherit = FALSE,
    type = "scatter3d",
    mode = "text"
  ) %>%
  layout(
    title = "Variance Surface of H² (Color = Mean H²) with Highlighted Path",
    scene = list(
      xaxis = list(title = "n_geno"),
      yaxis = list(title = "n_rep"),
      zaxis = list(title = "Var(H²)")
    )
  )

p







p <- plot_ly(
  x = x_vals,
  y = y_vals,
  z = var_mat,
  type = "surface",
  surfacecolor = mean_mat,
  colorscale = "Viridis",
  colorbar = list(title = "Mean H²")
)

# add the highlight line on the surface
p <- p %>% add_trace(
  data = highlight_pts,
  x = ~n_geno,
  y = ~n_rep,
  z = ~var_h2,
  type = "scatter3d",
  mode = "lines+markers",
  marker = list(size = 6, color = "red"),
  line   = list(color = "red", width = 5),
  name = "Selected Path"
)

# add mean(H2) text labels at each point
p <- p %>% add_text(
  data = highlight_pts,
  x = ~n_geno,
  y = ~n_rep,
  z = ~var_h2,
  text = ~sprintf("mean=%.3f", mean_h2),
  textposition = "top center",
  showlegend = FALSE
)

# layout
p <- p %>% layout(
  title = "Variance Surface of H² (Color = Mean H²) with Highlighted Path",
  scene = list(
    xaxis = list(title = "n_geno"),
    yaxis = list(title = "n_rep"),
    zaxis = list(title = "Var(H²)")
  )
)

p



##### 3D cloud plot #####
df_mean <- df_trait %>%
  group_by(n_geno, n_rep) %>%
  summarise(mean_h2 = mean(h2), .groups = "drop")

plot_ly(
  df_trait,
  x = ~n_geno,
  y = ~n_rep,
  z = ~h2,
  mode = "markers",
  marker = list(size = 3),
  type = "scatter3d"
) %>%
  layout(
    title = "3D Point Cloud of H² Distribution",
    scene = list(
      xaxis = list(title = "n_geno"),
      yaxis = list(title = "n_rep"),
      zaxis = list(title = "H² values")
    )
  )

plot_ly(df_trait,
        x = ~n_geno,
        y = ~n_rep,
        z = ~h2,
        mode = "markers",
        marker = list(size = 3, color = "lightgrey"),
        type = "scatter3d") %>%
  
  # add mean points
  add_markers(data = df_mean,
              x = ~n_geno,
              y = ~n_rep,
              z = ~mean_h2,
              marker = list(size = 6, color = "red"),
              name = "Mean H²") %>%
  
  layout(
    title = "3D Point Cloud of H² Distribution with Mean Values",
    scene = list(
      xaxis = list(title = "n_geno"),
      yaxis = list(title = "n_rep"),
      zaxis = list(title = "H² values")
    )
  )

#################  Stair Plot ########
df_trait <- df_trait %>%
  mutate(
    n_geno = factor(n_geno),
    n_rep = factor(n_rep)
  )

ggplot(df_trait, aes(x = n_rep, y = h2, group = n_rep)) +
  geom_violin(fill = "#4C9ED9", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.6) +   # optional
  facet_grid(n_geno ~ ., switch = "y") +  # stair-like layout
  labs(
    title = "Stair Plot of H² Distributions Across Genotype × Replicate Settings",
    x = "Number of Replicates",
    y = "H² Distribution",
    caption = "Each row = n_geno; Each violin = distribution of H² for a given n_rep"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "#e8e8e8"),
    strip.text.y.left = element_text(angle = 0),
    panel.spacing = unit(1.2, "lines"),
    plot.title = element_text(face = "bold")
  )

