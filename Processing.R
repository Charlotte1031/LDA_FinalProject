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
VisitData$SER_OS <- VisitData$S_OS + (VisitData$C_OS / 2)

# Right Eye (OD)
VisitData$S_OD <- rowMeans(VisitData[, c("AutoRef1SphOD", "AutoRef2SphOD", "AutoRef3SphOD")], na.rm = TRUE)
VisitData$C_OD <- rowMeans(VisitData[, c("AutoRef1CylOD", "AutoRef2CylOD", "AutoRef3CylOD")], na.rm = TRUE)
VisitData$SER_OD <- VisitData$S_OD + (VisitData$C_OD / 2)

# Step 2: Average SER across both eyes for each visit (i.e., time j)
VisitData$SER <- rowMeans(VisitData[, c("SER_OS", "SER_OD")], na.rm = TRUE)

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

# Step 5: Filter for patients with more than 5 visits without NA SER
FilteredData <- MergedData %>%
  filter(!is.na(SER)) %>%  # Remove rows where SER is NA
  group_by(PtID) %>%       # Group by PtID
  filter(n() > 5) %>%      # Keep only groups (patients) with more than 5 valid visits
  ungroup()                # Remove grouping to return to a regular data frame

# Get unique PtID for filtered dataset
UniquePtID <- unique(FilteredData$PtID)


# Part 2. Spagetti Plot

# Sample 10 unique patients from each treatment group
SampledPatients_trt <- FilteredData %>%
  filter(TrtGroup == "Atropine") %>%
  distinct(PtID, .keep_all = TRUE) %>%
  slice_sample(n = 10) %>%
  select(PtID)

SampledPatients_ctrl <- FilteredData %>%
  filter(TrtGroup == "Placebo") %>%
  distinct(PtID, .keep_all = TRUE) %>%
  slice_sample(n = 10) %>%
  select(PtID)

SampledPatients <- bind_rows(SampledPatients_trt, SampledPatients_ctrl)

# Filter the original data to include only records for the sampled patients
SampledData <- FilteredData %>%
  filter(PtID %in% SampledPatients$PtID)

# Set line types based on the treatment group
SampledData <- SampledData %>%
  mutate(LineType = ifelse(TrtGroup == "Atropine", "solid", "dashed"))

# Define the order of visits for consistent x-axis ordering
visit_order <- c("Enrollment", "Month 6 Visit", "Month 12 Visit", 
                 "Month 18 Visit", "Month 24 Visit", "Month 30 Visit")

# Convert Visit to a factor with specified order
SampledData <- SampledData %>%
  filter(Visit != "Run-in FU Randomization") %>%
  mutate(Visit = factor(Visit, levels = visit_order))

# Create the spaghetti plot with distinct line types and colors for each patient
ggplot(SampledData, aes(x = Visit, y = SER, group = PtID, color = as.factor(PtID))) +
  geom_line(aes(linetype = LineType), size = 0.5, alpha = 0.7) +  # Line for each patient with line type by group
  geom_point(size = 2, shape = 1) +                                        # Add points at each SER measurement
  labs(title = "SER Progression Over Time by Treatment Group",
       x = "Visit",
       y = "SER",
       color = "Patient ID",
       linetype = 'LineType') +
  scale_x_discrete(drop = FALSE) +                              # Ensure all visits appear on the x-axis
  scale_linetype_manual(values = c("solid", "dashed")) +        
  theme_minimal() +
  theme(legend.position = "none")                               # Optionally hide the legend if there are too many patients



# Calculate and add mean profiles
MeanProfileData <- MergedData %>%
  filter(Visit != "Run-in FU Randomization") %>%
  group_by(Visit, TrtGroup) %>%
  summarize(Mean_SER = mean(SER, na.rm = TRUE), .groups = 'drop')

# Combined plot with individual trajectories and mean profiles
ggplot() +
  # Individual trajectories
  geom_line(data = SampledData, 
            aes(x = Visit, y = SER, group = PtID, color = as.factor(PtID), 
                linetype = TrtGroup),
            size = 0.5, alpha = 0.7) +
  geom_point(data = SampledData,
             aes(x = Visit, y = SER, color = as.factor(PtID)),
             size = 2, shape = 1) +
  # Mean profiles
  geom_line(data = MeanProfileData,
            aes(x = Visit, y = Mean_SER, group = TrtGroup, linetype = TrtGroup),
            size = 1.2) +
  geom_point(data = MeanProfileData,
             aes(x = Visit, y = Mean_SER, color = TrtGroup),
             size = 3, shape = 21, fill = "white") +
  # Labels and theme
  labs(title = "SER Progression Over Time by Treatment Group",
       x = "Visit",
       y = "SER",
       color = "Patient ID",
       linetype = "Treatment Group") +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  theme(legend.position = "right")

### Displaying five patients

# only leave the last observation, which is the 5th visit
SampledData %>%
  select(
    PtID,
    SER,
    Sex,
    TreatmentGroup,
    Race,
    EyeColor,
    AgeAsofEnrollDt,
    is_mother_has_myopia,
    is_father_has_myopia,
    Visit
  ) %>%
  group_by(PtID) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  arrange(PtID) %>%
  select(-Visit) %>% # then rename the column , ageasenrolldt as age and Eyecolor as iris_color
  rename(age = AgeAsofEnrollDt, iris_color = EyeColor) %>%
  head(5) %>%
  knitr::kable()


# key summary statistics for demographics

MergedData %>% 
  select(PtID,
         Sex,
         TreatmentGroup,
         Race,
         EyeColor,
         AgeAsofEnrollDt,
         is_mother_has_myopia,
         is_father_has_myopia) %>%
  rename(age = AgeAsofEnrollDt, iris_color = EyeColor) %>%
  group_by(PtID) %>% 
  slice_tail(n = 1) %>%
  ungroup() %>%
  select(-PtID) -> Demographics

Demographics$age %>% summary()
Demographics %>%
  select(-age) %>%
  apply(.,2,table) %>%
  knitr::kable()

save(MergedData, file = "MergedData.RData")

#### clean dataset with selected variables

data <- MergedData[,c("PtID","Visit","SER","TrtGroup","MotherMyop","FatherMyop","Sex","Race","EyeColor","AgeAsofEnrollDt")]

mother_gene <- ifelse(data$MotherMyop=="Yes",1,0)
father_gene <- ifelse(data$FatherMyop=="Yes",1,0)
data$genetic <- mother_gene+father_gene

data$Age <- floor(data$AgeAsofEnrollDt)


# save as Rdata

save(data,file="data.Rdata")







