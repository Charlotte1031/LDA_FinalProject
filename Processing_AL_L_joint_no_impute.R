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

final_L_eye_data <- final_sep_eye_data %>% filter(Eye == "Left")

# Modify the clean data preparation to create AL_cat instead of groups
clean_data <- final_L_eye_data %>%
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
  # Create AL_cat based on baseline AL
  group_by(PtID) %>%
  mutate(
    Baseline_AL = first(AL_sp[Visit == "Month0"], default = NA),
    AL_cat = as.integer(Baseline_AL > 25)  # 0 for low AL (≤25), 1 for high AL (>25)
  ) %>%
  ungroup() %>%
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
  )

# Function to analyze data with joint model
analyze_joint_data <- function(clean_data) {
  # Relevel treatment group
  clean_data$TrtGroup <- relevel(as.factor(clean_data$TrtGroup), ref = "Placebo")
  
  # Fit joint model
  tryCatch({
    model <- lmer(AL_sp ~ Visit_numeric * as.factor(genetic) * TrtGroup * AL_cat +
                    Race + EyeColor + AgeAsofEnrollDt + 
                    (Visit_numeric | PtID) - Visit_numeric : as.factor(genetic) : TrtGroup : AL_cat,
                  data = clean_data,
                  na.action = na.omit)
    return(model)
  }, error = function(e) {
    print(paste("Error in model fitting:", e$message))
    return(NULL)
  })
}

# Function to extract model results
extract_joint_results <- function(model) {
  if(is.null(model)) {
    return(NULL)
  }
  
  # Extract fixed effects
  fixed_effects <- summary(model)$coefficients
  conf_int <- confint(model, method="Wald")
  
  # Create results dataframe
  results <- data.frame(
    term = rownames(fixed_effects),
    estimate = fixed_effects[,"Estimate"],
    std.error = fixed_effects[,"Std. Error"],
    statistic = fixed_effects[,"t value"]
  )
  
  # Add confidence intervals
  results$conf.low <- conf_int[rownames(fixed_effects), 1]
  results$conf.high <- conf_int[rownames(fixed_effects), 2]
  
  # Add p-values (using normal approximation)
  results$p.value <- 2 * (1 - pnorm(abs(results$statistic)))
  
  # Extract random effects
  random_effects <- as.data.frame(VarCorr(model))
  random_effects <- data.frame(
    group = random_effects$grp,
    var_name = random_effects$var1,
    variance = random_effects$vcov
  )
  
  return(list(
    fixed_effects = results,
    random_effects = random_effects
  ))
}

# Perform analysis
joint_model <- analyze_joint_data(clean_data)
results <- extract_joint_results(joint_model)

# Print results
if(!is.null(results)) {
  cat("\nFixed Effects:\n")
  print(results$fixed_effects, digits = 3)
  cat("\nRandom Effects:\n")
  print(results$random_effects, digits = 3)
}
