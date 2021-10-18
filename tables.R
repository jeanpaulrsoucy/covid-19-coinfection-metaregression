# Predictors and Microbiology of Co-Infection in Patients with COVID-19: Living Rapid Review Update and Meta-Regression #
# Langford et al. #
# Script author: Jean-Paul R. Soucy #

# Reproduce manuscript tables #

# Note: This script assumes the working directory is set to the root directory of the project
# This is most easily achieved by using the provided covid-19-coinfection-metaregression.Rproj in RStudio

# create output directory (if it doesn't exist already)
dir.create("tables", showWarnings = FALSE)

# load libraries

## meta-analysis and meta-regression
library(meta)
library(metafor)

## data processing
library(tibble)
library(dplyr)
library(broom)

# load functions
source("funs.R")

# load data
dat <- read.csv("data/coinfection.csv", header = TRUE, stringsAsFactors = FALSE)

# create variables

## bacterial infection variable (highest of coinfection and secondary infection)
dat <- dat %>%
  rowwise() %>%
  mutate(bacterial_infection = max(coinfection, secondary_infection, na.rm = TRUE)) %>%
  ungroup()

## factor variables
dat <- dat %>%
  mutate(
    setting = ifelse(Setting == "Outpatients", "Hospital/Outpatients", Setting),
    setting = factor(setting,
                     levels = c("Hospital/Outpatients", "Hospital", "Hospital ICU"),
                     labels = c("Hospital/Outpatient", "Hospital", "ICU")),
    body_site = factor(body_site, levels = c("Respiratory", "Blood", "Respiratory and Blood", "Unclear")),
    study_end_month = format(as.Date(End.Date, "%d/%m/%Y"), "%b"),
    study_end_month_grp = factor(
      case_when(
        study_end_month %in% month.abb[1:5] ~ "Jan–May",
        study_end_month %in% month.abb[6:10] ~ "Jun–Oct"
      ),
      levels = c("Jan–May", "Jun–Oct")
    ),
    region = factor(Region,
                    levels = c("North America", "Europe", "Asia",
                               "Middle East", "South and Central America", "Multiple")),
    risk_of_bias = factor(bias,
                          levels = c("low", "moderate", "high"),
                          labels = c("Low", "Moderate", "High"))
  )

## continuous variables
dat <- dat %>%
  mutate(
    age = as.numeric(Age),
    percent_female = as.integer(Female.patients..n.) / Sample * 100,
    percent_mechanical_vent = as.integer(Mechanical.Ventilation..n.) / Sample * 100,
    percent_ards = as.integer(ARDS..n.) / Sample * 100,
    percent_icu = as.integer(ICU..n.) / Sample * 100,
    percent_smoker = as.integer(Smokers..n.) / Sample * 100,
    percent_copd = as.integer(COPD..n.) / Sample * 100,
    percent_cvd = as.integer(Cardiovascular.disease..n.) / Sample * 100,
    percent_diabetes = as.integer(Diabetes..n.) / Sample * 100,
    percent_corticosteroid = as.integer(corticosteroids) / Sample * 100,
    percent_il6 = as.integer(il_6_inhibitors) / Sample * 100
  ) %>%
  rowwise() %>%
  mutate(percent_severe = ifelse(
    all(is.na(percent_mechanical_vent), is.na(percent_ards), is.na(percent_icu)),
    NA,
    max(percent_mechanical_vent, percent_ards, percent_icu, na.rm = TRUE))
    ) %>%
  ungroup()

# run models

## define outcomes
outcomes <- c("coinfection", "secondary_infection", "bacterial_infection")

## define adjustment variables
adjust <- c("age", "percent_severe")

## run models
vars <- c("setting", "study_end_month_grp", "region", "risk_of_bias",
          "age", "percent_female", "percent_mechanical_vent", "percent_smoker",
          "percent_copd", "percent_cvd", "percent_diabetes",
          "percent_corticosteroid", "percent_il6")

## create summary table
summary_table <- lapply(vars, function(v) {
  run_rma_mod(dat, outcomes, v, adjust) %>%
    unlist(recursive = FALSE) %>%
    extract_rma_mod(v) %>%
    make_table_row(v) %>%
    make_sample_columns(outcomes, v, adjust)
  }) %>%
  bind_rows() # combine all the rows

## replace some data (alternative adjustment sets)

## % Mechanical ventilation - don't control for percent_severe
summary_table <- replace_rows(dat, outcomes,
                              var = "percent_mechanical_vent",
                              adjust = "age",
                              summary_table = summary_table,
                              characteristic = "% Mechanical ventilation",
                              n_rows = 1)

## Setting - don't control for percent_severe
summary_table <- replace_rows(dat, outcomes,
                              var = "setting",
                              adjust = "age",
                              summary_table = summary_table,
                              characteristic = "Setting",
                              n_rows = 4)

## add table headers
summary_table <- summary_table %>%
  make_table_headers(outcomes)

## write summary table
write.table(summary_table,
            "tables/summary_table.csv", sep = ",",
            row.names = FALSE, col.names = FALSE)

# sensitivity analysis: unadjusted estimates using full sample size

## build table
summary_table_unadjusted_only <- lapply(vars, function(v) {
  run_rma_mod(dat, outcomes, v, adjust = NULL, adjusted = FALSE) %>%
    extract_rma_mod(v) %>%
    make_table_row(v) %>%
    make_sample_columns(outcomes, v, adjust = NULL, adjusted = FALSE)
}) %>%
  bind_rows() 

## add table headers
summary_table_unadjusted_only <- summary_table_unadjusted_only %>%
  make_table_headers(outcomes, adjusted = FALSE)

## write summary table
write.table(summary_table_unadjusted_only,
            "tables/summary_table_unadjusted_only.csv", sep = ",",
            row.names = FALSE, col.names = FALSE)

# sensitivity analysis: repeat metaregression table for co-infection
# excluding Type == "Not-specified"
## (in the primary analysis, we assume non-specified infections are co-infections)

## filter dataset
dat_exclude_not_specified <- dat %>%
  filter(Type != "Not-specified")

## create summary table
summary_table_exclude_not_specified <- lapply(vars, function(v) {
  run_rma_mod(dat_exclude_not_specified, outcomes, v, adjust) %>%
    unlist(recursive = FALSE) %>%
    extract_rma_mod(v) %>%
    make_table_row(v) %>%
    make_sample_columns(outcomes, v, adjust)
}) %>%
  bind_rows() # combine all the rows

## replace some data (alternative adjustment sets)

## % Mechanical ventilation - don't control for percent_severe
summary_table_exclude_not_specified <- replace_rows(dat_exclude_not_specified, outcomes,
                              var = "percent_mechanical_vent",
                              adjust = "age",
                              summary_table = summary_table_exclude_not_specified,
                              characteristic = "% Mechanical ventilation",
                              n_rows = 1)

## Setting - don't control for percent_severe
summary_table_exclude_not_specified <- replace_rows(dat_exclude_not_specified, outcomes,
                              var = "setting",
                              adjust = "age",
                              summary_table = summary_table_exclude_not_specified,
                              characteristic = "Setting",
                              n_rows = 4)

## add table headers
summary_table_exclude_not_specified <- summary_table_exclude_not_specified %>%
  make_table_headers(outcomes, exclude_not_specified = TRUE)

## write summary table
write.table(summary_table_exclude_not_specified,
            "tables/summary_table_exclude_not_specified.csv", sep = ",",
            row.names = FALSE, col.names = FALSE)
