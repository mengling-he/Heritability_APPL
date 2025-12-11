library(tidyr)
library(dplyr)
library(reshape2)
library(tibble)

library(lme4)
library(lmerTest)
library(ggplot2)

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))


# source the heritability calculation function
source("../../../Code/heritability_fun.R")


# ----read data----
features_df0_top_side <- read.csv("../data/side_top_features_thresholding_linearimputed_day.csv",stringsAsFactors = FALSE)
colnames(features_df0_top_side)
str(features_df0_top_side$Date)
var(features_df0_top_side[,c('side_hue_max')])

var(features_df_top_side[,c('side_hue_max')])
##### delete dead  #####
slow_growers <- strsplit(trimws(readLines("../data/slow_grower.txt")), ",")[[1]]
features_df_top_side<- features_df0_top_side  %>% filter(! Plant.Info %in% slow_growers)

##### generate a Day column  #####
features_df_top_side$Date <- as.Date(features_df_top_side$Date)
features_df_top_side$Day <- as.numeric(features_df_top_side$Date - min(features_df_top_side$Date))
colnames(features_df_top_side)

#------- model: time as fixed effect  -------
metric_cols <- setdiff(
  colnames(features_df_top_side),
  c("Plant.Info", "Date", "Genotype", "Replicate", "Day")
)
metric_cols_sub <- setdiff(
  metric_cols,
  c("side_hue_max"))

h2_results <- calc_H2_fixedtime(features_df_top_side,metric_cols_sub)
write.csv(h2_results,"../result/h2_timefixed_day_threshold.csv",row.names = FALSE)

h2_results_sort <- h2_results[order(-h2_results$H2), ]
sum(h2_results$H2 > 0.2)
highH2list <- subset(h2_results, H2 > 0.2)$Trait

ggplot(h2_results, aes(x = H2, y = Trait)) +
  geom_col(fill = "skyblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.y = element_text(size=2)) +
  labs(title = "H2 for different metric based on daily imputed data", x = "", y = "H2")

ggplot(h2_results_sort[1:sum(h2_results$H2 > 0.2),], 
       aes(x = reorder(Trait, H2), y = H2)) +
  geom_col(fill = "skyblue") +
  coord_flip() +  # makes it horizontal and easier to read
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 12)
  ) +
  labs(
    title = "Top H2 (> 0.2) based on daily imputed data",
    x = "Trait",
    y = "H2"
  )




#------- model: across time  -------
h2_results_timeseries <- calc_H2_timeseries(features_df_top_side,metric_cols_sub,time_col='Date')
h2_results_timeseries$Date <- as.Date(h2_results_timeseries$Time)
#write.csv(h2_results_timeseries,"../result/h2_timeseries_day_threshold.csv",row.names = FALSE)
filtered_H2_trait_timeseries <- h2_results_timeseries %>%
  group_by(Trait) %>%
  filter(any(H2 > 0.4)) %>%
  ungroup()
filtered_H2_trait_timeseries %>%
  summarise(n_unique_traits = n_distinct(Trait))

plot_heatmap(filtered_H2_trait_timeseries,columns_col = "Date",annot=FALSE,title = 'Heritability (H²>0.4) over Time')

ggplot(h2_results_timeseries %>%
         group_by(Trait) %>%
         filter(any(H2 > 0.5)) %>%
         ungroup(), aes(x = Date, y = H2, color = Trait, group = Trait)) +
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
    title = "Heritability (H²>0.5) over Time",
    x = "Date",
    y = expression(H^2),
    color = "Trait"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )








# ----use LDA without feature selection----
library(MASS)
# LDA formula: class ~ features
# delete features of top_blue_min and top_saturation_max
metric_cols_lda <- metric_cols_sub
df_scaled <- features_df_top_side
df_scaled[metric_cols_lda] <- scale(df_scaled[metric_cols_lda])

lda_model <- lda(Genotype ~ ., data = df_scaled[, c("Genotype", metric_cols_lda)])
lda_pred <- predict(lda_model, df_scaled)
lda_x <- data.frame(lda_pred$x) # LD1, LD2, ...

LDA_dim <- colnames(lda_x)[1:5]# select the first 5 dimensions
LDA_df <- cbind(df_scaled[c('Date','Genotype','Replicate')],lda_x[LDA_dim])
LDA_df$Date <- as.Date(LDA_df$Date)
LDA_df$Day <- as.numeric(LDA_df$Date - min(LDA_df$Date))

h2_results_lda <- calc_H2_fixedtime(LDA_df,LDA_dim,day_col = "Day")

ggplot(h2_results_lda,aes(x = reorder(Trait, H2), y = H2)) +
  geom_col(fill = "skyblue") +
  geom_text(aes(label = round(H2, 3)),   # display H² value rounded to 3 decimals
            hjust = 1.1,                # position text slightly outside the bar
            size = 3.5) +    
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.y = element_text(size=2)) +
  labs(title = "LDA H2 based on Daily Imputed Data - Threshold", x = "", y = "H2")

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

# 
# # ----use LDA without feature selection time series----
# h2_results_lda_timeseries <- calc_H2_timeseries(LDA_df,LDA_dim,time_col='Date')
# h2_results_lda_timeseries$Date <- as.Date(h2_results_lda_timeseries$Time)
# 
# plot_heatmap(h2_results_lda_timeseries,columns_col = "Date",annot=TRUE,title = 'Heritability  over Time')
# 
# ggplot(h2_results_lda_timeseries, aes(x = Date, y = H2, color = Trait, group = Trait)) +
#   geom_line(size = 1) +
#   geom_point(size = 2) +
#   geom_vline(
#     xintercept = as.Date('2024-06-19'),   # <-- adjust year as needed
#     linetype = "dashed",
#     color = "black",
#     linewidth = 0.8
#   ) +annotate(
#     "text",
#     x = as.Date("2024-06-19"),
#     y = 1.05,
#     label = "June 19",
#     angle = 0,
#     vjust = -0.5,
#     hjust = 1,
#     size = 3.5
#   ) +
#   labs(
#     title = "Heritability of LDA over Time",
#     x = "Date",
#     y = expression(H^2),
#     color = "Trait"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "right",
#     plot.title = element_text(face = "bold", hjust = 0.5)
#   )

