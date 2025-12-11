
# ----compare day and date group----
# 
# h2_day_imputed <- h2_results
# h2_group_imputed <- h2_results_fromimputed
# h2_combined <- merge(
#   h2_day_imputed[, c("Trait", "H2")],
#   h2_group_imputed[, c("Trait", "H2")],
#   by = "Trait",
#   suffixes = c("_day", "_group"),
#   all = TRUE
# )
# colnames(h2_combined) <- c("Trait","H2_daily_imputed","H2_group_imputed")
# 
# write.csv(
#   h2_combined,"results/H2_results/0_H2combine_day_group.csv",row.names = FALSE
# )
h2_combined <- read.csv("results/H2_results/0_H2combine_day_group.csv")
h2_combined_sub <- h2_combined %>%
  filter(
    !is.na(Trait),                    # remove NA trait row
    if_any(starts_with("H2"), ~ .x > 0.2)   # keep rows with ≥1 H2 > 0.4
  )
h2_combined_sub_long <- h2_combined_sub %>%
  pivot_longer(
    cols = -Trait,
    names_to = "Method",
    values_to = "H2"
  )
ggplot(h2_combined_sub_long, aes(x = Trait, y = H2, fill = Method)) +
  geom_col(position = position_dodge()) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Comparison of H² Across Methods",
       x = "Trait", y = "Heritability (H²)")

h2_combined_sub %>%
  filter(H2_daily_imputed > H2_group_imputed) %>%
  nrow()

# ----use LDA without feature selection----
library(MASS)
# LDA formula: class ~ features
# delete features of top_blue_min and top_saturation_max
lda_model <- lda(Genotype ~ ., data = df_scaled[, c("Genotype", metric_cols_sub2)])
lda_pred <- predict(lda_model, df_scaled)
lda_x <- data.frame(lda_pred$x) # LD1, LD2, ...

LDA_dim <- colnames(lda_x)[1:5]# select the first 10 dimensions
LDA_df <- cbind(df_scaled[c('Date','Genotype','Replicate')],lda_x[LDA_dim])
LDA_df$Date <- as.Date(LDA_df$Date)
h2_results_lda2 <- calc_H2_fixedtime(LDA_df,LDA_dim,day_col = "Date")

ggplot(h2_results_lda,aes(x = reorder(Trait, H2), y = H2)) +
  geom_col(fill = "skyblue") +
  geom_text(aes(label = round(H2, 3)),   # display H² value rounded to 3 decimals
            hjust = -0.1,                # position text slightly outside the bar
            size = 3.5) +    
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.y = element_text(size=2)) +
  labs(title = "LDA H2 based on Daily Imputed Data", x = "", y = "H2")

##### lda scaling#####################
lda_scale <- data.frame(lda_model$scaling)
lda_scale$Variable <- rownames(lda_scale)

# LD1 top 20
top_LD1 <- lda_scale %>%
  dplyr::select(Variable, LD1) %>%
  left_join(h2_results %>% dplyr::select(Trait, H2), by = c("Variable" = "Trait")) %>%
  arrange(desc(abs(LD1))) %>%
  slice(1:20)  # top 20 LD1 contributors

# LD2 top 20
top_LD2 <- lda_scale %>%
  dplyr::select(Variable, LD2) %>%
  left_join(h2_results %>% dplyr::select(Trait, H2), by = c("Variable" = "Trait")) %>%
  arrange(desc(abs(LD2))) %>%
  slice(1:20)  # top 20 LD2 contributors

ggplot(top_LD1, aes(x = reorder(Variable, LD1), y = LD1, fill = H2)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "plasma") +
  labs(
    title = "Top 20 Variable Contributions to LD1",
    x = "",
    y = "Contribution to LD1",
    fill = "H²"
  ) +
  theme_minimal()

ggplot(top_LD2, aes(x = reorder(Variable, LD2), y = LD2, fill = H2)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "plasma") +
  labs(
    title = "Top 20 Variable Contributions to LD2",
    x = "",
    y = "Contribution to LD2",
    fill = "H²"
  ) +
  theme_minimal()

# ggplot(top_LD2, aes(x = reorder(Variable, LD2), y = LD2)) +
#   geom_col(fill = "salmon") +
#   coord_flip() +
#   labs(title = "Top 20 Contributions to LD2", x = "", y = "Contribution to LD2") +
#   theme_minimal()
# ggplot(top_LD1, aes(x = reorder(Variable, LD1), y = LD1)) +
#   geom_col(fill = "skyblue") +
#   coord_flip() +
#   labs(title = "Top 20 Contributions to LD1", x = "", y = "Contribution to LD1") +
#   theme_minimal()

# ----use LDA without feature selection time series----
h2_results_lda_timeseries <- calc_H2_timeseries(LDA_df,LDA_dim,time_col='Date')
h2_results_lda_timeseries$Date <- as.Date(h2_results_lda_timeseries$Time)

plot_heatmap(h2_results_lda_timeseries,columns_col = "Date",annot=TRUE,title = 'Heritability  over Time')

ggplot(h2_results_lda_timeseries, aes(x = Date, y = H2, color = Trait, group = Trait)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_vline(
    xintercept = as.Date('2024-06-19'),   # <-- adjust year as needed
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +annotate(
    "text",
    x = as.Date("2024-06-19"),
    y = 1.05,
    label = "June 19",
    angle = 0,
    vjust = -0.5,
    hjust = 1,
    size = 3.5
  ) +
  labs(
    title = "Heritability of LDA over Time",
    x = "Date",
    y = expression(H^2),
    color = "Trait"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

