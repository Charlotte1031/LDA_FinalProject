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
      TRUE ~ "Other"
    ),
    Race = factor(Race, levels = c("Black", "White", "Other")),
    
    # Clean Visit names
    Visit = case_when(
      Visit == "Enrollment" ~ "Month0",
      Visit == "Run-in FU Randomization" ~ "Randomization",
      Visit == "Month 6 Visit" ~ "Month06",
      Visit == "Month 12 Visit" ~ "Month12", 
      Visit == "Month 18 Visit" ~ "Month18",
      Visit == "Month 24 Visit" ~ "Month24",
      Visit == "Month 30 Visit" ~ "Month30"
    )
  ) %>%
  # Assign baseline AL group based on the AL value at Month0
  group_by(PtID) %>%
  mutate(
    Baseline_AL = first(AL[Visit == "Month0"], default = NA),
    AL_group = cut(
      Baseline_AL,
      breaks = c(-Inf, 25, Inf),  # Breakpoints for groups
      labels = c("Group1", "Group2"),
      include.lowest = TRUE
    )
  ) %>%
  ungroup() %>%
  # Create binary indicator variables for each subgroup
  mutate(
    AL_group1 = as.integer(AL_group == "Group1"),
    AL_group2 = as.integer(AL_group == "Group2"),
    AL_group3 = as.integer(AL_group == "Group3"),
    AL_group4 = as.integer(AL_group == "Group4")
  ) %>%
  filter(Visit != "Randomization") %>% 
  mutate(
    Visit_numeric = case_when(
      Visit == "Month0" ~ 0,
      Visit == "Month06" ~ 6,
      Visit == "Month12" ~ 12,
      Visit == "Month18" ~ 18,
      Visit == "Month24" ~ 24,
      Visit == "Month30" ~ 30
    )
  ) %>%
  mutate(AL_per_visit = AL / (Visit_numeric+1))


# Convert to wide format for imputation
wide_data <- clean_data %>%
  select(PtID, Visit, AL, TrtGroup, MotherMyop, FatherMyop, Sex, Race, EyeColor, AgeAsofEnrollDt, genetic, Age) %>%
  pivot_wider(
    id_cols = c(PtID, TrtGroup, MotherMyop, FatherMyop, Sex, Race, EyeColor, AgeAsofEnrollDt, genetic, Age),
    names_from = Visit,
    values_from = AL,
    names_prefix = "AL_"
  )

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
pred_mat["AL_Month0", c("AL_Month06","AL_Month12","AL_Month18","AL_Month24","AL_Month30")] <- 0
pred_mat["AL_Month06", c("AL_Month12","AL_Month18","AL_Month24","AL_Month30")] <- 0
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


  long_data <- long_data %>%
    mutate(
      AL_group = cut(
        AL,
        breaks = c(-Inf, 25, Inf),
        labels = c("Group1", "Group2"),
        include.lowest = TRUE
      ),
      AL_group1 = as.integer(AL_group == "Group1"),
      AL_group2 = as.integer(AL_group == "Group2")
    )
  
  # Extract Group 1 and Group 2 datasets
  group1_data <- long_data %>%
    filter(AL_group == "Group1")
  
  group2_data <- long_data %>%
    filter(AL_group == "Group2")
  
  # Fit models for each group
  models <- list()
  
  # Try-catch blocks to handle potential convergence issues
  tryCatch({
    models$group1 <- lmer(AL ~ Visit * as.factor(genetic) * TrtGroup + 
                            Race + EyeColor + AgeAsofEnrollDt +
                            (1 | PtID),
                          data = group1_data,
                          na.action = na.omit)
  }, error = function(e) {
    print(paste("Error in Group 1, imputation", i, ":", e$message))
    return(NULL)
  })
  
  tryCatch({
    models$group2 <- lmer(AL ~ Visit * as.factor(genetic) * TrtGroup + 
                            Race + EyeColor + AgeAsofEnrollDt +
                            (1 | PtID),
                          data = group2_data,
                          na.action = na.omit)
  }, error = function(e) {
    print(paste("Error in Group 2, imputation", i, ":", e$message))
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
        Visit == "Month06" ~ 6,
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

save(all_models, file = "all_models.RData")

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

# Combine all results
all_results <- rbind(results_group1, results_group2)

# Print results
print(all_results, digits = 3)
save(all_results,file="result.Rdata")

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

####### model diagnostic
# Residual Plot for a Single Model
residuals_plot <- function(model, title = "Residuals vs Fitted") {
  residuals <- residuals(model)
  fitted <- fitted(model)
  
  plot(fitted, residuals, 
       main = title, 
       xlab = "Fitted Values", 
       ylab = "Residuals", 
       pch = 16, col = "blue")
  abline(h = 0, col = "red", lty = 2)
}

# Example for group1 model in the first imputation
residuals_plot(all_models[[1]]$group1, title = "Residuals vs Fitted: Group 1")


# QQ Plot for Residuals
qq_plot_residuals <- function(model, title = "QQ Plot for Residuals") {
  qqnorm(residuals(model), main = title, pch = 16, col = "blue")
  qqline(residuals(model), col = "red", lty = 2)
}

# Example for group1 model
qq_plot_residuals(all_models[[1]]$group1, title = "QQ Plot: Group 1 Residuals")

# QQ Plot for Random Effects
qq_plot_random <- function(model, title = "QQ Plot for Random Effects") {
  ranef_data <- ranef(model)$PtID[, "(Intercept)"]
  qqnorm(ranef_data, main = title, pch = 16, col = "blue")
  qqline(ranef_data, col = "red", lty = 2)
}

# Example for group1 model
qq_plot_random(all_models[[1]]$group1, title = "QQ Plot: Group 1 Random Effects")

# install.packages("performance")
# library(performance)
# 
# # Check diagnostics for group1 model
# check_model(all_models[[1]]$group1,check = c("qq", "linearity","reqq"))
# check_model(all_models[[1]]$group2,check = c("qq", "linearity","reqq"))
# check_model(all_models[[1]]$group3,check = c("qq", "linearity","reqq"))

# library(car)
# # VIF for group1 model
# vif(all_models[[1]]$group3)

# Extract random effects for the intercept
random_effects <- ranef(all_models[[1]]$group1)$PtID[, "(Intercept)"]
shapiro.test(random_effects)

# Residuals vs Fitted Plot for each group
par(mfrow = c(2, 1)) 

# Group 1 - 105
if (!is.null(group1_model)) {
  residuals_plot(group1_model, title = "Residuals vs Fitted: Group 1")
}

# Group 2 - 127
if (!is.null(group2_model)) {
  residuals_plot(group2_model, title = "Residuals vs Fitted: Group 2")
}


# QQ Plot for Residuals for each group
# Group 1
if (!is.null(group1_model)) {
  qq_plot_residuals(group1_model, title = "QQ Plot: Group 1 Residuals")
}

# Group 2
if (!is.null(group2_model)) {
  qq_plot_residuals(group2_model, title = "QQ Plot: Group 2 Residuals")
}

# QQ Plot for Random Effects for each group
# Group 1
if (!is.null(group1_model)) {
  qq_plot_random(group1_model, title = "QQ Plot: Group 1 Random Effects")
}

# Group 2
if (!is.null(group2_model)) {
  qq_plot_random(group2_model, title = "QQ Plot: Group 2 Random Effects")
}


dev.off()

