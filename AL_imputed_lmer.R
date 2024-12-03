load("data.RData")
library(tidyverse)
library(lme4)
library(nlme)
library(mice)


# Clean data as before, but focus on AL
clean_data <- data %>%
  # First reclassify Race
  mutate(
    Race = case_when(
      Race == "Black/African American" ~ "Black",
      Race == "White" ~ "White",
      TRUE ~ "Other"  # All other categories combined
    ),
    # Convert Race to factor with specific reference level
    Race = factor(Race, levels = c("Black", "White", "Other")),
    # Then clean Visit
    Visit = case_when(
      Visit == "Enrollment" ~ "Month0",
      Visit == "Run-in FU Randomization" ~ "Randomization",
      Visit == "Month 6 Visit" ~ "Month6",
      Visit == "Month 12 Visit" ~ "Month12", 
      Visit == "Month 18 Visit" ~ "Month18",
      Visit == "Month 24 Visit" ~ "Month24",
      Visit == "Month 30 Visit" ~ "Month30"
    )
  )

# Convert to wide format for imputation
wide_data <- clean_data %>%
  select(PtID, Visit, AL, TrtGroup, MotherMyop, FatherMyop, Sex, Race, EyeColor, AgeAsofEnrollDt, genetic, Age) %>%
  pivot_wider(
    id_cols = c(PtID, TrtGroup, MotherMyop, FatherMyop, Sex, Race, EyeColor, AgeAsofEnrollDt, genetic, Age),
    names_from = Visit,
    values_from = AL,
    names_prefix = "AL_"
  ) %>% 
  select(-AL_Randomization)

# Check missingness
missing_pattern <- wide_data %>% 
  select(starts_with("AL_")) %>% 
  apply(2, function(x) sum(is.na(x)))

print("Missing values by timepoint:")
print(missing_pattern)

# Set up imputation
library(mice)

# Create prediction matrix
pred_mat <- make.predictorMatrix(wide_data)

# Modify prediction matrix for temporal ordering
pred_mat["AL_Month0", c("AL_Month6","AL_Month12","AL_Month18","AL_Month24","AL_Month30")] <- 0
pred_mat["AL_Month6", c("AL_Month12","AL_Month18","AL_Month24","AL_Month30")] <- 0
pred_mat["AL_Month12", c("AL_Month18","AL_Month24","AL_Month30")] <- 0
pred_mat["AL_Month18", c("AL_Month24","AL_Month30")] <- 0
pred_mat["AL_Month24", "AL_Month30"] <- 0

# Set imputation methods
meth <- make.method(wide_data)
meth[grep("AL_", names(meth))] <- "pmm"  # Use predictive mean matching for AL

# Perform multiple imputation
imp <- mice(wide_data, 
            m = 100,           # Generate 20 imputed datasets
            maxit = 20,       # Number of iterations
            method = meth,    
            predictorMatrix = pred_mat,
            seed = 123)

# Check imputation quality
densityplot(imp)
plot(imp)





# Function to analyze one imputed dataset
analyze_imputed_dataset <- function(imp_data, i) {
  print(paste("Analyzing imputation", i))
  
  # Convert to long format
  long_data <- to_long(imp_data)
  
  # Create AL groups and indicator variables
  long_data <- long_data %>%
    mutate(
      AL_group = cut(
        AL,
        breaks = c(-Inf, 24.5, 26, 27.5, Inf),
        labels = c("Group1", "Group2", "Group3", "Group4"),
        include.lowest = TRUE
      ),
      AL_group1 = as.integer(AL_group == "Group1"),
      AL_group2 = as.integer(AL_group == "Group2"),
      AL_group3 = as.integer(AL_group == "Group3"),
      AL_group4 = as.integer(AL_group == "Group4")
    )
  
  # Split into groups
  group1_data <- long_data %>% 
    filter(AL_group1 == 1) %>% 
    select(-AL_group1, -AL_group2, -AL_group3, -AL_group4)
  
  group2_data <- long_data %>% 
    filter(AL_group2 == 1) %>% 
    select(-AL_group1, -AL_group2, -AL_group3, -AL_group4)
  
  group3_data <- long_data %>% 
    filter(AL_group3 == 1) %>% 
    select(-AL_group1, -AL_group2, -AL_group3, -AL_group4)
  
  # Fit models for each group
  models <- list()
  
  # Try-catch blocks to handle potential convergence issues
  tryCatch({
    models$group1 <- lmer(AL ~ Visit_numeric * as.factor(genetic) * TrtGroup + 
                            Race + EyeColor + AgeAsofEnrollDt +
                            (Visit_numeric | PtID),
                          data = group1_data,
                          na.action = na.omit)
  }, error = function(e) {
    print(paste("Error in Group 1, imputation", i, ":", e$message))
    return(NULL)
  })
  
  tryCatch({
    models$group2 <- lmer(AL ~ Visit_numeric * as.factor(genetic) * TrtGroup + 
                            Race + EyeColor + AgeAsofEnrollDt +
                            (Visit_numeric | PtID),
                          data = group2_data,
                          na.action = na.omit)
  }, error = function(e) {
    print(paste("Error in Group 2, imputation", i, ":", e$message))
    return(NULL)
  })
  
  tryCatch({
    models$group3 <- lmer(AL ~ Visit_numeric * as.factor(genetic) * TrtGroup + 
                            Race+ EyeColor + AgeAsofEnrollDt +
                            (Visit_numeric | PtID),
                          data = group3_data,
                          na.action = na.omit)
  }, error = function(e) {
    print(paste("Error in Group 3, imputation", i, ":", e$message))
    return(NULL)
  })
  
  return(models)
}

to_long <- function(data) {
  data %>%
    pivot_longer(
      cols = starts_with("AL_"),
      names_to = "Visit",
      values_to = "AL",
      names_prefix = "AL_"
    ) %>%
    mutate(
      Visit_numeric = case_when(
        Visit == "Month0" ~ 0,
        Visit == "Month6" ~ 6,
        Visit == "Month12" ~ 12,
        Visit == "Month18" ~ 18,
        Visit == "Month24" ~ 24,
        Visit == "Month30" ~ 30
      )
    )
}

get_coef_se <- function(model) {
  coef <- fixef(model)
  vcov <- as.matrix(vcov(model))
  se <- sqrt(diag(vcov))
  # Create numeric matrix instead of data frame
  matrix(c(coef, se), ncol = 2, 
         dimnames = list(names(coef), c("estimate", "std.error")))
}

# Analyze all imputed datasets
all_models <- list()
for(i in 1:100) {
  imp_data <- complete(imp, i)
  all_models[[i]] <- analyze_imputed_dataset(imp_data, i)
}

# Function to extract coefficients and SE for each group
extract_group_results <- function(models_list, group_name) {
  # Extract models for specific group
  group_models <- lapply(models_list, function(x) x[[group_name]])
  
  # Remove NULL models (if any failed to converge)
  group_models <- group_models[!sapply(group_models, is.null)]
  
  # Extract coefficients and standard errors
  results_list <- lapply(group_models, get_coef_se)
  
  # Combine results for Rubin's rules
  estimates <- do.call(rbind, lapply(results_list, function(x) x[,"estimate"]))
  std.errors <- do.call(rbind, lapply(results_list, function(x) x[,"std.error"]))
  
  # Apply Rubin's rules
  n_imputations <- nrow(estimates)
  pooled_estimates <- colMeans(estimates)
  within_var <- colMeans(std.errors^2)
  between_var <- apply(estimates, 2, var)
  total_var <- within_var + (1 + 1/n_imputations) * between_var
  pooled_se <- sqrt(total_var)
  
  # Create results dataframe
  results <- data.frame(
    term = colnames(estimates),
    estimate = pooled_estimates,
    std.error = pooled_se,
    statistic = pooled_estimates/pooled_se,
    p.value = 2 * (1 - pnorm(abs(pooled_estimates/pooled_se))),
    conf.low = pooled_estimates - 1.96 * pooled_se,
    conf.high = pooled_estimates + 1.96 * pooled_se,
    group = group_name
  )
  
  return(results)
}

# Get results for each group
results_group1 <- extract_group_results(all_models, "group1")
results_group2 <- extract_group_results(all_models, "group2")
results_group3 <- extract_group_results(all_models, "group3")

# Combine all results
all_results <- rbind(results_group1, results_group2, results_group3)

# Print results
print(all_results, digits = 3)

# Visualize key coefficients across groups
library(ggplot2)

ggplot(all_results %>%
         filter(term %in% c("TrtGroup", "Visit_numeric", "genetic", 
                            "Visit_numeric:genetic")) %>%
         mutate(term = factor(term, 
                              levels = c("TrtGroup", "Visit_numeric", 
                                         "genetic", "Visit_numeric:genetic"))), 
       aes(x = group, y = estimate, color = term)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = 0.5),
                width = 0.2) +
  theme_minimal() +
  coord_flip() +
  labs(title = "Key Coefficient Estimates by AL Group",
       x = "AL Group",
       y = "Estimate with 95% CI",
       color = "Parameter")


