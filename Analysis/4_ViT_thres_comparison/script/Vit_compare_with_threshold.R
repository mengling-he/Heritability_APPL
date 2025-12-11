library(tidyr)
library(dplyr)
library(naniar)
library(ggplot2)
library(lme4)

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))


# source the heritability calculation function
source("../../../Code/heritability_fun.R")



#----------Compare fixed time results--------

H2_fixed_threshold <- read.csv("../data/h2_timefixed_day_threshold.csv",stringsAsFactors = FALSE)
H2_fixed_vit<- read.csv("../data/h2_timefixed_day_vit.csv",stringsAsFactors = FALSE)


H2_fixed_vit$Dataset       <- "ViT"
H2_fixed_threshold$Dataset <- "Threshold"

h2_fixed_all <- rbind(H2_fixed_vit, H2_fixed_threshold)

# Make factors with nice order
h2_fixed_all$Dataset <- factor(h2_fixed_all$Dataset, levels = c("ViT", "Threshold"))

ggplot(h2_fixed_all,
       aes(x = Trait, y = H2, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Comparison of H²",
    x     = "Trait",
    y     = "H²"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )


######   plot the subset of comparison ######
traits_keep <- h2_fixed_all %>%
  group_by(Trait) %>%
  summarise(max_H2 = max(H2, na.rm = TRUE), .groups = "drop") %>%
  filter(max_H2 > 0.2) %>%
  pull(Trait)


h2_filtered <- h2_fixed_all %>%
  filter(Trait %in% traits_keep) %>%
  mutate(
    Trait = factor(
      Trait,
      levels = traits_keep[order(
        sapply(traits_keep, function(t)
          max(h2_fixed_all$H2[h2_fixed_all$Trait == t], na.rm = TRUE))
      )]
    )
  )
ggplot(h2_filtered,
       aes(x = reorder(Trait, H2), y = H2, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Comparison of H² (traits with max H² > 0.2 in at least one dataset)",
    x     = "Trait",
    y     = "H²"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1),
    panel.grid.minor = element_blank()
  )


###### get the traits that are higher in Vit#########
traits_vit_higher <- h2_fixed_all %>%
  select(Trait, Dataset, H2) %>%
  pivot_wider(names_from = Dataset, values_from = H2) %>%
  filter(ViT > Threshold) %>%
  pull(Trait)
h2_vit_higher <- h2_fixed_all %>%
  filter(Trait %in% traits_vit_higher)
ggplot(h2_vit_higher,
       aes(x = Trait, y = H2, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8)) +
  theme_minimal() +  
  labs(
    title = "ViT > Threshold",
    x     = "Traits",
    y     = "H²"
  ) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1))

###### get the traits that are higher in threshold#########
traits_threshold_higher <- h2_fixed_all %>%
  select(Trait, Dataset, H2) %>%
  pivot_wider(names_from = Dataset, values_from = H2) %>%
  filter(Threshold > ViT) %>%
  pull(Trait)
h2_threshold_higher <- h2_fixed_all %>%
  filter(Trait %in% traits_threshold_higher)
ggplot(h2_threshold_higher,
       aes(x = Trait, y = H2, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(
    title = "Threshold > ViT",
    x     = "Traits",
    y     = "H²"
  ) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1))









#----------compare LDA results --------

##### Read data and delet no variance traits#########
slow_growers <- strsplit(trimws(readLines("../data/slow_grower.txt")), ",")[[1]]

meta_cols = c("Plant.Info","Date", "Genotype", "Replicate")

# delete the 0-variance columns
df_vit <- read.csv("../data/side_top_features_vit_linearimputed_day.csv",stringsAsFactors = FALSE)
df_vit<- df_vit  %>% filter(! Plant.Info %in% slow_growers)
metric_all_vit  <- setdiff(colnames(df_vit), meta_cols)
var_vec <- sapply(df_vit[,metric_all_vit], var, na.rm = TRUE)
zero_var_cols <- names(var_vec[var_vec == 0])
zero_var_cols
df_vit_clean <- df_vit[, !(colnames(df_vit) %in% zero_var_cols)]



df_threshold <- read.csv("../data/side_top_features_thresholding_linearimputed_day.csv",stringsAsFactors = FALSE)
df_threshold<- df_threshold  %>% filter(! Plant.Info %in% slow_growers)
metric_all_threshold  <- setdiff(colnames(df_threshold), meta_cols)
var_vec <- sapply(df_threshold[,metric_all_threshold], var, na.rm = TRUE)
zero_var_cols <- names(var_vec[var_vec == 0])
zero_var_cols
df_threshold_clean <- df_threshold[, !(colnames(df_threshold) %in% zero_var_cols)]



##### Pipeline-------
run_lda_sets <- function(df,
                         meta_cols = c("Plant.Info", "Date", "Genotype", "Replicate"),
                         n_dim = 5,
                         day_col = "Day") {
  
  metric_all  <- setdiff(colnames(df), meta_cols)
  metric_top  <- grep("^top_",  metric_all, value = TRUE)
  metric_side <- grep("^side_", metric_all, value = TRUE)
  
  metric_sets <- list(
    all  = metric_all,
    top  = metric_top,
    side = metric_side
  )
  
  ## store h2 results for all three sets
  out_h2 <- lapply(names(metric_sets), function(set_name) {
    traits <- metric_sets[[set_name]]
    
    if (length(traits) == 0) return(NULL)
    
    message("Running LDA for set: ", set_name, " (", length(traits), " traits)")
    
    ## ---- scale ----
    X_scaled <- as.data.frame(scale(df[, traits, drop = FALSE]))
    colnames(X_scaled) <- traits
    
    ## ---- rebuild data ----
    df_scaled <- cbind(df[, meta_cols, drop = FALSE], X_scaled)
    df_scaled$Date <- as.Date(df_scaled$Date)
    
    ## ---- run LDA ----
    formula_lda <- as.formula(
      paste("Genotype ~", paste(traits, collapse = " + "))
    )
    
    lda_model <- lda(formula_lda, data = df_scaled)
    lda_pred  <- predict(lda_model, df_scaled)
    lda_x     <- as.data.frame(lda_pred$x)
    
    ## ---- first k LDs ----
    keep_dims <- head(colnames(lda_x), n_dim)
    
    ## ---- final LDA dataframe ----
    LDA_df <- cbind(
      df_scaled[, meta_cols, drop = FALSE],
      lda_x[, keep_dims, drop = FALSE]
    )
    
    LDA_df[[day_col]] <- as.numeric(LDA_df$Date - min(LDA_df$Date, na.rm = TRUE))
    
    ## ---- calculate H2 ----
    h2_df <- calc_H2_fixedtime(LDA_df, keep_dims, day_col = day_col)
    
    ## add the Set name
    h2_df$Set <- set_name
    
    return(h2_df)
  })
  
  ## ---- combine into one dataframe ----
  out_h2 <- do.call(rbind, out_h2)
  rownames(out_h2) <- NULL
  
  return(out_h2)
}


##### run pipeline --------
# run pipeline on each
res_vit <- run_lda_sets(df_vit_clean)
res_threshold <- run_lda_sets(df_threshold_clean)
dim(res_vit)
dim(res_threshold)


##### Analyze result --------
res_vit$Dataset       <- "ViT"
res_threshold$Dataset <- "Threshold"

h2_all <- rbind(res_vit, res_threshold)

# Make factors with nice order
h2_all$Trait   <- factor(h2_all$Trait, levels = paste0("LD", 1:5))
h2_all$Set     <- factor(h2_all$Set, levels = c("all", "top", "side"))
h2_all$Dataset <- factor(h2_all$Dataset, levels = c("ViT", "Threshold"))

h2_all$Group <- interaction(h2_all$Dataset, h2_all$Set, sep = "_")

ggplot(h2_all,
       aes(x = Group, y = H2, fill = Set)) +
  geom_col() +
  facet_wrap(~ Trait, nrow = 1) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "H² by LDA dimension",
    x     = "Dataset",
    y     = "H²"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ggplot(h2_all,
#        aes(x = Trait, y = H2, fill = Set)) +
#   geom_col(position = position_dodge(width = 0.8)) +
#   facet_wrap(~ Dataset) +
#   scale_fill_brewer(palette = "Set2") +
#   labs(
#     title = "Comparison of H² for LDA dimensions",
#     x     = "LDA Dimension",
#     y     = "H²"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid.minor = element_blank()
#   )

##### LD1 partition---------
# keep only LD1
ld1_df <- h2_all %>%
  filter(Trait == "LD1") 

ld1_long <- ld1_df %>%
  mutate(
    Total_Var = Var_Genotype + Var_Replicate + Var_Residual
  ) %>%
  select(Group, Var_Genotype, Var_Replicate, Var_Residual, Total_Var) %>%
  pivot_longer(
    cols = c(Var_Genotype, Var_Replicate, Var_Residual),
    names_to = "Component",
    values_to = "Variance"
  ) %>%
  mutate(
    Component = dplyr::recode(
      Component,
      Var_Genotype  = "Genotype",
      Var_Replicate = "Genotype:Replicate",
      Var_Residual  = "Residual"
    ),
    Prop      = Variance / Total_Var,
    pct_label = scales::percent(Prop, accuracy = 0.1)  # 比如 "43%"
  )

ggplot(ld1_long, aes(x = Group, y = Variance, fill = Component)) +
  geom_col() +
  geom_text(
    aes(label = pct_label),
    position = position_stack(vjust = 0.5),  
    color = "white",
    size = 3
  ) +
  labs(
    title = "Variance decomposition for LD1",
    x     = "Dataset",
    y     = "Variance",
    fill  = "Component"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )
