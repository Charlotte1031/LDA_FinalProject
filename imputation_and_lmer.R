library(tidyverse)
library(lme4)
library(broom.mixed)

load("data.RData")

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
    ),
    log_SER = log(-SER)
  )

# 2. 转换为宽格式
wide_data <- clean_data %>%
  select(PtID, Visit,log_SER, SER, TrtGroup, MotherMyop, FatherMyop, Sex, Race, EyeColor, AgeAsofEnrollDt, genetic, Age) %>%
  pivot_wider(
    id_cols = c(PtID, TrtGroup, MotherMyop, FatherMyop, Sex, Race, EyeColor, AgeAsofEnrollDt, genetic, Age),
    names_from = Visit,
    values_from = SER,
    names_prefix = "SER_"
  ) %>% select(-SER_Randomization)

view(wide_data)

# check missingness

wide_data %>% select(
  SER_Month6,
  SER_Month12,
  SER_Month18,
  SER_Month24,
  SER_Month30
) %>% apply(.,2,function(x) sum(is.na(x)))


# start imputing missing values

library(mice)

# install naniar if not installed
if (!requireNamespace("naniar", quietly = TRUE)) {
  install.packages("naniar")
}
library(naniar)

# visualize missing data

md.pattern(wide_data)

vis_miss(wide_data)


# Create prediction matrix for longitudinal data

# 1. Set up prediction matrix
pred_mat <- make.predictorMatrix(wide_data)

# 2. Modify prediction matrix to ensure temporal ordering
# Set future time points to 0 to prevent them from being used in prediction
# This ensures we only use past and present data to impute missing values
pred_mat["SER_Month0", c("SER_Month6","SER_Month12","SER_Month18","SER_Month24","SER_Month30")] <- 0
pred_mat["SER_Month6", c("SER_Month12","SER_Month18","SER_Month24","SER_Month30")] <- 0
pred_mat["SER_Month12", c("SER_Month18","SER_Month24","SER_Month30")] <- 0
pred_mat["SER_Month18", c("SER_Month24","SER_Month30")] <- 0
pred_mat["SER_Month24", "SER_Month30"] <- 0

# 3. Set imputation methods
meth <- make.method(wide_data)
# Set method for all SER columns to normal distribution
# Since SER is a continuous variable, we use normal distribution for imputation
meth[grep("SER_", names(meth))] <- "pmm"  

# 4. Perform multiple imputation
imp <- mice(wide_data, 
            m = 20,           # Generate 10 imputed datasets
            maxit = 50,      # Number of iterations
            method = meth,   # Specified imputation methods
            predictorMatrix = pred_mat,  # Custom prediction matrix
            seed = 123)      # Set seed for reproducibility

# 5. Check imputation quality
# Examine distributions of imputed values versus observed values
densityplot(imp)
# Check convergence of imputation algorithm
plot(imp)

# 6. Get completed datasets
# Get all imputed datasets in long format
completed_data <- complete(imp, action = "long")  
# Or get a single imputed dataset
completed_data_1 <- complete(imp, 1)  # Get first imputed dataset


view(completed_data_1)

# Starting Lmer model for imputed data and use Rubin Fromula
to_long <- function(data) {
  data %>%
    pivot_longer(
      cols = starts_with("SER_"),
      names_to = "Visit",
      values_to = "SER",
      names_prefix = "SER_"
    ) %>%
    # Convert Visit to actual month numbers
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

# 2. Create list to store models
model_list <- vector("list", 10)  

# 3. Loop through each imputation
for(i in 1:20) {
  imp_data <- complete(imp, i)
  long_data <- to_long(imp_data)
  
  model_list[[i]] <- lmer(SER ~ TrtGroup +genetic + Visit + Sex + Race + EyeColor +
                            AgeAsofEnrollDt + (1 | PtID),
                          data = long_data,
                          REML = TRUE)
}

get_coef_se <- function(model) {
  coef <- fixef(model)
  vcov <- as.matrix(vcov(model))
  se <- sqrt(diag(vcov))
  # Create numeric matrix instead of data frame
  matrix(c(coef, se), ncol = 2, 
         dimnames = list(names(coef), c("estimate", "std.error")))
}

# 2. Extract results from all models
results_list <- lapply(model_list, get_coef_se)

# 3. Combine results for pooling
estimates <- do.call(rbind, lapply(results_list, function(x) x[,"estimate"]))
std.errors <- do.call(rbind, lapply(results_list, function(x) x[,"std.error"]))

# 4. Pool results using Rubin's rules manually
n_imputations <- length(model_list)
pooled_estimates <- colMeans(estimates)
within_var <- colMeans(std.errors^2)
between_var <- apply(estimates, 2, var)
total_var <- within_var + (1 + 1/n_imputations) * between_var
pooled_se <- sqrt(total_var)

# 5. Create summary table
final_results <- data.frame(
  term = colnames(estimates),
  estimate = pooled_estimates,
  std.error = pooled_se,
  statistic = pooled_estimates/pooled_se,
  p.value = 2 * (1 - pnorm(abs(pooled_estimates/pooled_se))),
  conf.low = pooled_estimates - 1.96 * pooled_se,
  conf.high = pooled_estimates + 1.96 * pooled_se
)

# 6. Print results in nice format
print(final_results, digits = 3)





