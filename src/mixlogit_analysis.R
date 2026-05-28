# ============================================================================
# Mixed Logit Model Analysis for BWS Conjoint Data
# ============================================================================
# This script performs mixed logit model estimation on BWS survey data
# Output: Overall coefficients/WTP and individual-level coefficients/WTP

# ============================================================================
# 1. Load Libraries
# ============================================================================

library(mlogit)
library(tidyverse)
library(dplyr)
library(readr)

# ============================================================================
# 2. Data Loading and Preparation
# ============================================================================

# Load test data
df <- read.csv("data/test_data.csv")

# Display basic info
cat("Data dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
cat("Respondents:", length(unique(df$ID)), "\n")
cat("Tasks per respondent:", max(df$task), "\n")
cat("Column names:\n")
print(names(df))

# ============================================================================
# 3. Data Transformation for mlogit
# ============================================================================

# mlogit requires a specific data format
# Create choice ID: respondent_task_type combination
df$choice_id <- paste(df$ID, df$task, df$type, sep="_")

# Prepare attributes for the model
# Define attribute columns (exclude ID, cluster_id, task, prof_id, alt, choice, type, choice_id)
attribute_cols <- c(
  "project_forest",
  "project_sea_grass",
  "cobenefit_promotion",
  "cobenefit_safety",
  "fores_promo",
  "fores_safe",
  "sea_promo",
  "sea_safe",
  "price"
)

# Create mlogit format data
# mlogit requires: choice, id, alternative, and attributes
mlogit_data <- df %>%
  select(ID, choice_id, alt, choice, all_of(attribute_cols)) %>%
  arrange(choice_id, alt)

# Convert to mlogit format
mlogit_data <- mlogit.data(
  mlogit_data,
  choice = "choice",
  shape = "long",
  alt.var = "alt",
  id.var = c("ID", "choice_id")
)

cat("\nmlogit data prepared.\n")
cat("Observations:", nrow(mlogit_data), "\n")
cat("Choices:", nrow(mlogit_data) / 4, "\n")  # 4 alternatives per choice

# ============================================================================
# 4. Mixed Logit Model Estimation
# ============================================================================

cat("\n=== ESTIMATING MIXED LOGIT MODEL ===\n")
cat("Random parameters: project_forest, project_sea_grass, cobenefit_promotion,\n")
cat("                   cobenefit_safety, fores_promo, fores_safe, sea_promo, sea_safe\n")
cat("Fixed parameters: price, research (not included in current model)\n")

# Model formula
# Note: research is not included as it's only for experimental design
formula_mixlogit <- choice ~ 
  project_forest + 
  project_sea_grass + 
  cobenefit_promotion + 
  cobenefit_safety + 
  fores_promo + 
  fores_safe + 
  sea_promo + 
  sea_safe + 
  price

# Estimate mixed logit model
# Random parameters: normal distribution (except price - fixed)
mixlogit_model <- mlogit(
  formula_mixlogit,
  data = mlogit_data,
  random = list(
    project_forest = 1,
    project_sea_grass = 1,
    cobenefit_promotion = 1,
    cobenefit_safety = 1,
    fores_promo = 1,
    fores_safe = 1,
    sea_promo = 1,
    sea_safe = 1
  ),
  rpar = c(
    project_forest = "n",
    project_sea_grass = "n",
    cobenefit_promotion = "n",
    cobenefit_safety = "n",
    fores_promo = "n",
    fores_safe = "n",
    sea_promo = "n",
    sea_safe = "n"
  ),
  panel = TRUE,
  halton = NA,
  draws = 100
)

# Display model summary
cat("\n=== MODEL SUMMARY ===\n")
print(summary(mixlogit_model))

# ============================================================================
# 5. Extract Overall Coefficients
# ============================================================================

cat("\n=== EXTRACTING OVERALL COEFFICIENTS ===\n")

# Get coefficient summary
coef_summary <- summary(mixlogit_model)$CoefTable
overall_coefs <- as.data.frame(coef_summary)
overall_coefs$Parameter <- rownames(overall_coefs)
rownames(overall_coefs) <- NULL

# Rename columns for clarity
colnames(overall_coefs) <- c("Estimate", "StdError", "tStat", "pValue", "Parameter")

# Separate mean and SD
overall_coefs_mean <- overall_coefs %>%
  filter(!grepl("sd\\.", Parameter)) %>%
  select(Parameter, Estimate, StdError, tStat, pValue) %>%
  mutate(Type = "Mean")

overall_coefs_sd <- overall_coefs %>%
  filter(grepl("sd\\.", Parameter)) %>%
  select(Parameter, Estimate) %>%
  mutate(Type = "SD") %>%
  rename(Parameter_SD = Parameter)

# Clean parameter names
overall_coefs_mean$Parameter <- gsub("^.*:", "", overall_coefs_mean$Parameter)

cat("Overall Coefficients (Mean) and Standard Deviations:\n")
print(overall_coefs)

# ============================================================================
# 6. Calculate Overall WTP
# ============================================================================

cat("\n=== CALCULATING OVERALL WTP ===\n")

# Extract price coefficient (negative sign for WTP calculation)
price_coef <- coef(mixlogit_model)["price"]
cat("Price coefficient:", price_coef, "\n")

# Extract mean estimates for random parameters
random_params <- c(
  "project_forest", 
  "project_sea_grass", 
  "cobenefit_promotion",
  "cobenefit_safety", 
  "fores_promo", 
  "fores_safe", 
  "sea_promo", 
  "sea_safe"
)

# WTP = -coefficient / price_coefficient
# (negative because we want positive WTP values for positive attributes)
wtp_overall <- data.frame(
  Parameter = random_params,
  Coefficient = coef(mixlogit_model)[paste0(random_params, ".mean")],
  PriceCoefficient = price_coef,
  stringsAsFactors = FALSE
)

wtp_overall$WTP <- -wtp_overall$Coefficient / wtp_overall$PriceCoefficient

# Add fixed parameters
wtp_overall$Parameter_Clean <- wtp_overall$Parameter
wtp_overall <- wtp_overall %>%
  select(Parameter_Clean, Coefficient, WTP) %>%
  arrange(desc(WTP))

colnames(wtp_overall) <- c("Attribute", "Coefficient", "WTP")

cat("\nOverall WTP (Willingness to Pay):\n")
cat("(Unit: % price increase for each unit increase in attribute)\n")
print(wtp_overall)

# ============================================================================
# 7. Calculate Individual-Level Coefficients
# ============================================================================

cat("\n=== CALCULATING INDIVIDUAL-LEVEL COEFFICIENTS ===\n")

# Get individual-level coefficients using conditional expectation
# This requires extracting the predicted random coefficients for each respondent

# Predict individual-level coefficients
individual_coefs <- rpar(mixlogit_model, newdata = mlogit_data, drawing = TRUE)

cat("Individual coefficient predictions generated.\n")
cat("Dimensions:", nrow(individual_coefs), "x", ncol(individual_coefs), "\n")

# Aggregate by respondent (average across choices if needed)
individual_summary <- as.data.frame(individual_coefs)
individual_summary$ID <- row.names(individual_summary)
individual_summary <- individual_summary %>%
  select(ID, everything()) %>%
  mutate(ID = as.numeric(ID))

# If there are duplicate respondents, take mean
individual_summary <- individual_summary %>%
  group_by(ID) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup()

cat("\nIndividual-level coefficients summary (first 10 respondents):\n")
print(head(individual_summary, 10))

# ============================================================================
# 8. Calculate Individual-Level WTP
# ============================================================================

cat("\n=== CALCULATING INDIVIDUAL-LEVEL WTP ===\n")

# WTP for each individual = -coefficient / price_coefficient
individual_wtp <- individual_summary %>%
  mutate(
    project_forest_wtp = -project_forest / price,
    project_sea_grass_wtp = -project_sea_grass / price,
    cobenefit_promotion_wtp = -cobenefit_promotion / price,
    cobenefit_safety_wtp = -cobenefit_safety / price,
    fores_promo_wtp = -fores_promo / price,
    fores_safe_wtp = -fores_safe / price,
    sea_promo_wtp = -sea_promo / price,
    sea_safe_wtp = -sea_safe / price
  )

cat("\nIndividual-level WTP summary (first 10 respondents):\n")
print(head(individual_wtp, 10))

# ============================================================================
# 9. Save Results to CSV Files
# ============================================================================

cat("\n=== SAVING RESULTS ===\n")

# Create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Save overall coefficients
write.csv(
  overall_coefs,
  "results/overall_coefficients.csv",
  row.names = FALSE
)
cat("✓ Saved: results/overall_coefficients.csv\n")

# Save overall WTP
write.csv(
  wtp_overall,
  "results/overall_wtp.csv",
  row.names = FALSE
)
cat("✓ Saved: results/overall_wtp.csv\n")

# Save individual coefficients
write.csv(
  individual_summary,
  "results/individual_coefficients.csv",
  row.names = FALSE
)
cat("✓ Saved: results/individual_coefficients.csv\n")

# Save individual WTP
write.csv(
  individual_wtp,
  "results/individual_wtp.csv",
  row.names = FALSE
)
cat("✓ Saved: results/individual_wtp.csv\n")

# ============================================================================
# 10. Summary Statistics and Visualization
# ============================================================================

cat("\n=== SUMMARY STATISTICS FOR INDIVIDUAL WTP ===\n")

wtp_attributes <- individual_wtp %>%
  select(ends_with("_wtp")) %>%
  rename_all(~gsub("_wtp", "", .))

summary_stats <- data.frame(
  Attribute = names(wtp_attributes),
  Mean = colMeans(wtp_attributes),
  SD = apply(wtp_attributes, 2, sd),
  Min = apply(wtp_attributes, 2, min),
  Max = apply(wtp_attributes, 2, max),
  Median = apply(wtp_attributes, 2, median)
)

print(summary_stats)

# Save summary statistics
write.csv(
  summary_stats,
  "results/wtp_summary_statistics.csv",
  row.names = FALSE
)
cat("\n✓ Saved: results/wtp_summary_statistics.csv\n")

# ============================================================================
# 11. Create Visualization (Optional)
# ============================================================================

cat("\n=== CREATING VISUALIZATIONS ===\n")

# Box plot of individual WTP by attribute
png("results/wtp_distribution_boxplot.png", width = 1000, height = 600)
boxplot(wtp_attributes, 
        main = "Distribution of Individual-Level WTP by Attribute",
        ylab = "WTP (% price increase)",
        xlab = "Attributes",
        las = 2)
dev.off()
cat("✓ Saved: results/wtp_distribution_boxplot.png\n")

# Bar plot of mean WTP
png("results/wtp_mean_barplot.png", width = 800, height = 600)
barplot(summary_stats$Mean, 
        names.arg = summary_stats$Attribute,
        main = "Mean WTP by Attribute",
        ylab = "Mean WTP (% price increase)",
        las = 2,
        col = "steelblue")
dev.off()
cat("✓ Saved: results/wtp_mean_barplot.png\n")

# ============================================================================
# 12. Final Summary
# ============================================================================

cat("\n")
cat("=" %*% 70, "\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("=" %*% 70, "\n")
cat("\nOutput files saved in 'results/' directory:\n")
cat("  1. overall_coefficients.csv - Overall model coefficients\n")
cat("  2. overall_wtp.csv - Overall WTP estimates\n")
cat("  3. individual_coefficients.csv - Individual-level coefficients (100 respondents)\n")
cat("  4. individual_wtp.csv - Individual-level WTP estimates (100 respondents)\n")
cat("  5. wtp_summary_statistics.csv - Summary statistics of individual WTP\n")
cat("  6. wtp_distribution_boxplot.png - Visualization of WTP distribution\n")
cat("  7. wtp_mean_barplot.png - Mean WTP by attribute\n")
cat("\n")
cat("Model Specification:\n")
cat("  - Random parameters (8): project_forest, project_sea_grass, cobenefit_promotion,\n")
cat("                           cobenefit_safety, fores_promo, fores_safe, sea_promo, sea_safe\n")
cat("  - Fixed parameters (1): price\n")
cat("  - Distribution: Normal for all random parameters\n")
cat("  - Sample size: ", length(unique(df$ID)), " respondents\n")
cat("  - Observations: ", nrow(df), " choice situations\n")
cat("\n")
