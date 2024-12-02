library(tidyverse)
library(lme4)
library(broom.mixed)

load("data.RData")



# 1. First clean and reshape the data
clean_data <- data %>%
  # Reclassify Race
  mutate(
    Race = case_when(
      Race == "Black/African American" ~ "Black",
      Race == "White" ~ "White",
      TRUE ~ "Other"  
    ),
    Race = factor(Race, levels = c("Black", "White", "Other")),
    # Clean Visit
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

# 2. Create baseline dataset
baseline_data <- clean_data %>%
  filter(Visit == "Month0") %>%
  select(PtID, SER) %>%
  rename(baseline_SER = SER)

# 3. Calculate changes from baseline
long_data <- clean_data %>%
  filter(Visit != "Month0", Visit != "Randomization") %>%  # Remove baseline and randomization visits
  left_join(baseline_data, by = "PtID") %>%
  mutate(SER_change = SER - baseline_SER)  # Calculate change from baseline

# 4. Convert to wide format for imputation
wide_data <- long_data %>%
  select(PtID, Visit, SER_change, TrtGroup, MotherMyop, FatherMyop, Sex, Race, 
         EyeColor, AgeAsofEnrollDt, genetic, Age) %>%
  pivot_wider(
    id_cols = c(PtID, TrtGroup, MotherMyop, FatherMyop, Sex, Race, 
                EyeColor, AgeAsofEnrollDt, genetic, Age),
    names_from = Visit,
    values_from = SER_change,
    names_prefix = "SER_change_"
  )

# 5. Set up prediction matrix for imputation
pred_mat <- make.predictorMatrix(wide_data)

# Modify prediction matrix for temporal ordering
pred_mat["SER_change_Month6", c("SER_change_Month12","SER_change_Month18",
                                "SER_change_Month24","SER_change_Month30")] <- 0
pred_mat["SER_change_Month12", c("SER_change_Month18","SER_change_Month24",
                                 "SER_change_Month30")] <- 0
pred_mat["SER_change_Month18", c("SER_change_Month24","SER_change_Month30")] <- 0
pred_mat["SER_change_Month24", "SER_change_Month30"] <- 0

# 6. Perform multiple imputation
meth <- make.method(wide_data)
meth[grep("SER_change_", names(meth))] <- "pmm"

imp <- mice(wide_data, 
            m = 20,
            maxit = 50,
            method = meth,
            predictorMatrix = pred_mat,
            seed = 123)

# 7. Convert back to long format for modeling
to_long <- function(data) {
  data %>%
    pivot_longer(
      cols = starts_with("SER_change_"),
      names_to = "Visit",
      values_to = "SER_change",
      names_prefix = "SER_change_"
    ) %>%
    mutate(
      Visit_numeric = case_when(
        Visit == "Month6" ~ 6,
        Visit == "Month12" ~ 12,
        Visit == "Month18" ~ 18,
        Visit == "Month24" ~ 24,
        Visit == "Month30" ~ 30
      )
    )
}

# 8. Fit models
model_list <- vector("list", 20)
for(i in 1:20) {
  imp_data <- complete(imp, i)
  long_data <- to_long(imp_data)
  
  model_list[[i]] <- lmer(SER_change ~ TrtGroup * genetic + Visit + 
                            Sex + Race + EyeColor + AgeAsofEnrollDt + (1 | PtID),
                          data = long_data,
                          REML = TRUE)
}

# 9. Pool results
get_coef_se <- function(model) {
  coef <- fixef(model)
  vcov <- as.matrix(vcov(model))
  se <- sqrt(diag(vcov))
  matrix(c(coef, se), ncol = 2, 
         dimnames = list(names(coef), c("estimate", "std.error")))
}

results_list <- lapply(model_list, get_coef_se)
estimates <- do.call(rbind, lapply(results_list, function(x) x[,"estimate"]))
std.errors <- do.call(rbind, lapply(results_list, function(x) x[,"std.error"]))

n_imputations <- length(model_list)
pooled_estimates <- colMeans(estimates)
within_var <- colMeans(std.errors^2)
between_var <- apply(estimates, 2, var)
total_var <- within_var + (1 + 1/n_imputations) * between_var
pooled_se <- sqrt(total_var)

final_results <- data.frame(
  term = colnames(estimates),
  estimate = pooled_estimates,
  std.error = pooled_se,
  statistic = pooled_estimates/pooled_se,
  p.value = 2 * (1 - pnorm(abs(pooled_estimates/pooled_se))),
  conf.low = pooled_estimates - 1.96 * pooled_se,
  conf.high = pooled_estimates + 1.96 * pooled_se
)

print(final_results, digits = 3)