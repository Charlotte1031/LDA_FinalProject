library(dplyr)
library(ggplot2)

# Part A. Data Import and merge

# Load VisitData and PtRoster datasets
VisitData <- read.table("Data Tables/MTS1ClinicTesting.txt", 
                        header = TRUE, sep = "|", fill = TRUE, stringsAsFactors = FALSE)
head(VisitData)

PtRoster <- read.table("Data Tables/MTS1PtRoster.txt", 
                       header = TRUE, sep = "|", fill = TRUE, stringsAsFactors = FALSE)
head(PtRoster)

# Load the additional dataset for parental myopia and other demographic information
DemoMedOcuHx <- read.table("Data Tables/MTS1DemoMedOcuHx.txt", 
                           header = TRUE, sep = "|", fill = TRUE, stringsAsFactors = FALSE)
head(DemoMedOcuHx)

#TODO: Handling baseline and run-in randomization 
# VisitData$Visit == 'Run-in FU Randomization'
# VisitData$Visit == 'Enrollment'
# if enrollment 

# Step 1: Calculate SER for each eye
# Left Eye (OS)
VisitData$S_OS <- rowMeans(VisitData[, c("AutoRef1SphOS", "AutoRef2SphOS", "AutoRef3SphOS")], na.rm = TRUE)
VisitData$C_OS <- rowMeans(VisitData[, c("AutoRef1CylOS", "AutoRef2CylOS", "AutoRef3CylOS")], na.rm = TRUE)
VisitData$A_OS <- rowMeans(VisitData[, c("AxialLen1OS", "AxialLen2OS", "AxialLen3OS")], na.rm = TRUE)

VisitData$SER_OS <- VisitData$S_OS + (VisitData$C_OS / 2)

# Right Eye (OD)
VisitData$S_OD <- rowMeans(VisitData[, c("AutoRef1SphOD", "AutoRef2SphOD", "AutoRef3SphOD")], na.rm = TRUE)
VisitData$C_OD <- rowMeans(VisitData[, c("AutoRef1CylOD", "AutoRef2CylOD", "AutoRef3CylOD")], na.rm = TRUE)
VisitData$A_OD <- rowMeans(VisitData[, c("AxialLen1OD", "AxialLen2OD", "AxialLen3OD")], na.rm = TRUE)

VisitData$SER_OD <- VisitData$S_OD + (VisitData$C_OD / 2)

# Step 2: Average SER across both eyes for each visit (i.e., time j)
VisitData$SER <- rowMeans(VisitData[, c("SER_OS", "SER_OD")], na.rm = TRUE)
VisitData$AL <- rowMeans(VisitData[, c("A_OD", "A_OS")], na.rm = TRUE)
# Step 3: Merge data with PtRoster and DemoMedOcuHx datasets
# Merging VisitData with PtRoster on Patient ID
MergedData <- merge(VisitData, PtRoster, by = "PtID")

# Merging the result with DemoMedOcuHx on Patient ID
MergedData <- merge(MergedData, DemoMedOcuHx, by = "PtID")

# Step 4: Create necessary variables
# Time (Visit) is already in VisitData as 'Visit'
# Group/Treatment assignment
MergedData$TreatmentGroup <- MergedData$TrtGroup

# Parental myopia status (investigating variables)
MergedData$is_mother_has_myopia <- MergedData$MotherMyop
MergedData$is_father_has_myopia <- MergedData$Father

library(tidyr)
# Prepare separate eye dataset
sep_eye_data <- MergedData %>%
  pivot_longer(
    cols = c(SER_OS, SER_OD),
    names_to = "Eye",
    values_to = "SER_sp"
  ) %>%
  pivot_longer(
    cols = c(A_OS, A_OD),
    names_to = "AL_Eye",
    values_to = "AL_sp"
  ) %>%
  mutate(
    # Create eye-specific PtID
    PtID = case_when(
      Eye == "SER_OS" ~ paste0(PtID, "_L"),
      Eye == "SER_OD" ~ paste0(PtID, "_R")
    ),
    # Ensure eye matches for SER and AL
    Eye = case_when(
      Eye == "SER_OS" ~ "Left",
      Eye == "SER_OD" ~ "Right"
    )
  ) %>%
  filter(
    (Eye == "Left" & AL_Eye == "A_OS") | 
      (Eye == "Right" & AL_Eye == "A_OD")
  ) %>%
  select(-AL_Eye)

# Create genetic variable
sep_eye_data$mother_gene <- ifelse(sep_eye_data$MotherMyop=="Yes", 1, 0)
sep_eye_data$father_gene <- ifelse(sep_eye_data$FatherMyop=="Yes", 1, 0)
sep_eye_data$genetic <- sep_eye_data$mother_gene + sep_eye_data$father_gene

# Age calculation
sep_eye_data$Age <- floor(sep_eye_data$AgeAsofEnrollDt)

# Select final columns similar to original data processing
final_sep_eye_data <- sep_eye_data %>%
  select(
    PtID, 
    Visit, 
    SER, 
    SER_sp,
    AL,
    AL_sp,
    TrtGroup, 
    MotherMyop, 
    FatherMyop, 
    Sex, 
    Race, 
    EyeColor, 
    AgeAsofEnrollDt,
    Eye,
    genetic,
    Age
  )

# Save the dataset
save(final_sep_eye_data, file = "final_sep_eye_data.RData")

# Load required libraries
library(tidyverse)
library(mice)
library(lme4)

# Load the separated eye-level data
load("final_sep_eye_data.RData")

# Clean data as before, but focus on AL
clean_data <- final_sep_eye_data %>%
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
      Visit == "Month 6 Visit" ~ "Month6",
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
      breaks = c(-Inf, 24.5, 25, 25.5, Inf),  # Breakpoints for groups
      labels = c("Group1", "Group2", "Group3", "Group4"),
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
      Visit == "Month6" ~ 6,
      Visit == "Month12" ~ 12,
      Visit == "Month18" ~ 18,
      Visit == "Month24" ~ 24,
      Visit == "Month30" ~ 30
    )
  ) %>%
  mutate(AL_per_visit = AL / (Visit_numeric+1))

# Convert to wide format for imputation
wide_data <- clean_data %>%
  select(PtID, Visit, AL_sp, TrtGroup, MotherMyop, FatherMyop, Sex, Race, EyeColor, AgeAsofEnrollDt, genetic, Age, Eye) %>%
  pivot_wider(
    id_cols = c(PtID, TrtGroup, MotherMyop, FatherMyop, Sex, Race, EyeColor, AgeAsofEnrollDt, genetic, Age, Eye),
    names_from = Visit,
    values_from = AL_sp,
    names_prefix = "AL_"
  ) 

# Check missingness
missing_pattern <- wide_data %>% 
  select(starts_with("AL_")) %>% 
  apply(2, function(x) sum(is.na(x)))

print("Missing values by timepoint:")
print(missing_pattern)

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
            m = 100,           # Generate 100 imputed datasets
            maxit = 20,        # Number of iterations
            method = meth,    
            predictorMatrix = pred_mat,
            seed = 123)

# Check imputation quality
densityplot(imp)
plot(imp)

# Function to convert back to long format
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

# Analyze imputed datasets
analyze_imputed_dataset <- function(imp_data, i) {
  print(paste("Analyzing imputation", i))
  
  # Convert to long format
  long_data <- to_long(imp_data)
  
  # Create AL groups and indicator variables
  long_data <- long_data %>%
    mutate(
      AL_group = cut(
        AL,
        breaks = c(-Inf, 24.5, 25, 25.5, Inf),
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
  
  group4_data <- long_data %>% 
    filter(AL_group4 == 1) %>% 
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
                            Race + EyeColor + AgeAsofEnrollDt +
                            (Visit_numeric | PtID),
                          data = group3_data,
                          na.action = na.omit)
  }, error = function(e) {
    print(paste("Error in Group 3, imputation", i, ":", e$message))
    return(NULL)
  })
  
  tryCatch({
    models$group4 <- lmer(AL ~ Visit_numeric * as.factor(genetic) * TrtGroup + 
                            Race + EyeColor + AgeAsofEnrollDt +
                            (Visit_numeric | PtID),
                          data = group3_data,
                          na.action = na.omit)
  }, error = function(e) {
    print(paste("Error in Group 4, imputation", i, ":", e$message))
    return(NULL)
  })
  
  return(models)
}

# Function to extract coefficients and standard errors
get_coef_se <- function(model) {
  coef <- fixef(model)
  vcov <- as.matrix(vcov(model))
  se <- sqrt(diag(vcov))
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
results_group4 <- extract_group_results(all_models, "group4")

# Combine all results
all_results <- rbind(results_group1, results_group2, results_group3, results_group4)

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
