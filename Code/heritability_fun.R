
library(lme4)
library(dplyr)
library(tidyr)

library(ggplot2)
library(reshape2)
library(RColorBrewer)


###########################calc_H2_1dim#########################################
calc_H2 <- function(data, traits) {
  # the model is use Genotype and Genotype:Replicate as random effect
  results <- list()
  
  for (trait in traits) {
    # Fit model
    model <- lmer(
      as.formula(paste(trait, "~ (1|Genotype)")),
      data = data
    )
    
    # Extract variance components
    var_comp <- as.data.frame(VarCorr(model)) %>%
      dplyr::select(grp, vcov) %>%
      tidyr::pivot_wider(names_from = grp, values_from = vcov)
    
    # Ensure missing components become 0
    if (!"Genotype" %in% colnames(var_comp)) var_comp$Genotype <- 0
    if (!"Residual" %in% colnames(var_comp)) var_comp$Residual <- 0
    
    # Variance components
    sigma_G2 <- var_comp$Genotype
    sigma_E2 <- var_comp$Residual
    
    # Broad-sense H²
    h2 <- round(sigma_G2 / (sigma_G2 + sigma_E2), 4)
    
    # Store results
    results[[trait]] <- data.frame(
      Trait = trait,
      Var_Genotype = sigma_G2,
      Var_Residual = var_comp$Residual,
      H2 = h2
    )
  }
  
  # Combine all results into one dataframe
  h2_results <- do.call(rbind, results)
  rownames(h2_results) <- NULL
  
  return(h2_results)
}

###########################calc_H2_fixedtime#########################################
calc_H2_fixedtime <- function(data, traits, day_col = "Day") {
# the model is use Genotype and Genotype:Replicate as random effect
  results <- list()
  
  for (trait in traits) {
    message("Processing: ", trait)
    
    # Fit model
    model <- lmer(
      as.formula(paste(trait, "~", day_col, "+ (1|Genotype) + (1|Genotype:Replicate)")),
      data = data
    )
    
    # Extract variance components
    var_comp <- as.data.frame(VarCorr(model)) %>%
      dplyr::select(grp, vcov) %>%
      tidyr::pivot_wider(names_from = grp, values_from = vcov)
    
    # Ensure missing components become 0
    if (!"Genotype" %in% colnames(var_comp)) var_comp$Genotype <- 0
    if (!"Genotype:Replicate" %in% colnames(var_comp)) var_comp$`Genotype:Replicate` <- 0
    if (!"Residual" %in% colnames(var_comp)) var_comp$Residual <- 0
    
    # Variance components
    sigma_G2 <- var_comp$Genotype
    sigma_E2 <- var_comp$`Genotype:Replicate` + var_comp$Residual
    
    # Broad-sense H²
    h2 <- round(sigma_G2 / (sigma_G2 + sigma_E2), 4)
    
    # Store results
    results[[trait]] <- data.frame(
      Trait = trait,
      Var_Genotype = sigma_G2,
      Var_Replicate = var_comp$`Genotype:Replicate`,
      Var_Residual = var_comp$Residual,
      H2 = h2
    )
  }
  
  # Combine all results into one dataframe
  h2_results <- do.call(rbind, results)
  rownames(h2_results) <- NULL
  
  return(h2_results)
}

###########################calc_H2_timeseries#########################################
calc_H2_timeseries <- function(data, traits, time_col = "Day", time_groups = NULL) {
  if (is.null(time_groups)) {
    time_groups <- sort(unique(data[[time_col]]))
  }
  
  results_timeseries <- list()
  i <- 0
  
  for (d in time_groups) {
    i <- i + 1
    message("Processing Day: ", d)
    
    df_day <- subset(data, data[[time_col]] == d)
    day_results <- list()
    
    for (trait in traits) {
      cat("  Trait:", trait, "\n")
      
      # Fit LMM: (1|Genotype)
      model <- tryCatch(
        lmer(as.formula(paste(trait, "~ (1|Genotype)")), data = df_day),
        error = function(e) {
          warning("Model failed for trait ", trait, " on day ", d)
          return(NULL)
        }
      )
      
      if (is.null(model)) next
      
      # Extract variance components
      var_comp <- as.data.frame(VarCorr(model)) %>%
        dplyr::select(grp, vcov) %>%
        tidyr::pivot_wider(names_from = grp, values_from = vcov)
      
      # Ensure missing variance terms are 0
      if (!"Genotype" %in% names(var_comp)) var_comp$Genotype <- 0
      if (!"Residual" %in% names(var_comp)) var_comp$Residual <- 0
      
      # Variance components
      sigma_G2 <- var_comp$Genotype
      sigma_E2 <- var_comp$Residual
      
      # Broad-sense heritability
      h2 <- round(sigma_G2 / (sigma_G2 + sigma_E2), 4)
      
      # Store result for this trait and day
      day_results[[trait]] <- data.frame(
        Time = d,
        Trait = trait,
        Var_Genotype = sigma_G2,
        Var_Residual = sigma_E2,
        H2 = h2
      )
    }
    
    # Combine results for this date
    results_timeseries[[i]] <- do.call(rbind, day_results)
  }
  
  # Combine all days
  h2_df_timeseries <- do.call(rbind, results_timeseries)
  rownames(h2_df_timeseries) <- NULL
  
  return(h2_df_timeseries)
}

