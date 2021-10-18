# Predictors and Microbiology of Co-Infection in Patients with COVID-19: Living Rapid Review Update and Meta-Regression #
# Langford et al. #
# Script author: Jean-Paul R. Soucy #

# Reproduce manuscript figures #

# Note: This script assumes the working directory is set to the root directory of the project
# This is most easily achieved by using the provided covid-19-coinfection-metaregression.Rproj in RStudio

# create output directory (if it doesn't exist already)
dir.create("figures", showWarnings = FALSE)

# load libraries

## meta-analysis and meta-regression
library(meta)
library(metafor)

## data processing
library(dplyr)
library(lubridate)

## plotting
library(ggplot2)

# load functions
source("funs.R")

# load data
dat <- read.csv("data/coinfection.csv", header = TRUE, stringsAsFactors = FALSE)

# process data

## create bacterial infection variable (highest of coinfection and secondary infection)
dat <- dat %>%
  rowwise() %>%
  mutate(bacterial_infection = max(coinfection, secondary_infection, na.rm = TRUE)) %>%
  ungroup()

## modify "Setting" data (merge "Outpatients" (n = 1) with "Hospital/Outpatients")
dat$Setting <- ifelse(dat$Setting == "Outpatients", "Hospital/Outpatients", dat$Setting)

## create/modify stratifying variables for primary forest plots
dat <- dat %>%
  mutate(
    Setting = factor(Setting, levels = c("Hospital/Outpatients", "Hospital", "Hospital ICU"), labels = c("Hospital/Outpatient", "Hospital", "ICU")),
    `Body site` = factor(body_site, levels = c("Respiratory", "Blood", "Respiratory and Blood", "Unclear")),
    `Clinical definition` = factor(clinical_def, levels = c("Y", "N"), labels = c("Yes", "No"))
  )

## extract study end month
dat <- dat %>%
  mutate(
    end_month_int = month(as.Date(End.Date, "%d/%m/%Y")),
    end_year_int = year(as.Date(End.Date, "%d/%m/%Y")), # unnecessary since all end dates are in 2020
    end_date = as.Date(paste(end_year_int, end_month_int, "01", sep = "-"), "%Y-%m-%d"),
    `End month` = ifelse(is.na(end_date), "Not specified", format(end_date, "%b")),
    `End month` = factor(`End month`, levels = c(month.abb, "Not specified")),
    `End month` = factor(
      case_when(
        `End month` %in% c(month.abb[1:3]) ~ "Jan to Mar",
        `End month` %in% c(month.abb[4:6])  ~ "Apr to Jun",
        `End month` %in% c(month.abb[7:9]) ~ "Jul to Sep",
        `End month` %in% c(month.abb[10:12]) ~ "Oct to Dec",
        `End month` %in% c("Not specified") ~ "Not specified"
      ),
      levels = c("Jan to Mar", "Apr to Jun", "Jul to Sep", "Oct to Dec", "Not specified"))
  )

# primary forest plots (coinfection, secondary infection, bacterial infection)

## coinfection

### subset where coinfection has a value
dat_co <- dat %>%
  filter(!is.na(coinfection))

### setting
forest_plot(forest_calc(dat_co, "coinfection", "Setting"), xmax = 100,
            out_pdf = "figures/forest_coinfection_setting.pdf",
            out_width = 10, out_height = 5)

### body site
forest_plot(forest_calc(dat_co, "coinfection", "Body site"), xmax = 100,
            out_pdf = "figures/forest_coinfection_bodysite.pdf",
            out_width = 10, out_height = 5)

## clinical definition
forest_plot(forest_calc(dat_co, "coinfection", "Clinical definition"), xmax = 100,
            out_pdf = "figures/forest_coinfection_clinicaldef.pdf",
            out_width = 10, out_height = 4)

## secondary infection

### subset where secondary infection has a value
dat_si <- dat %>%
  filter(!is.na(secondary_infection))

### setting
forest_plot(forest_calc(dat_si, "secondary_infection", "Setting"), xmax = 100,
            out_pdf = "figures/forest_secondaryinfection_setting.pdf",
            out_width = 10, out_height = 5)

### body site
forest_plot(forest_calc(dat_si, "secondary_infection", "Body site"), xmax = 100,
            out_pdf = "figures/forest_secondaryinfection_bodysite.pdf",
            out_width = 10, out_height = 5)

## clinical definition
forest_plot(forest_calc(dat_si, "secondary_infection", "Clinical definition"), xmax = 100,
            out_pdf = "figures/forest_secondaryinfection_clinicaldef.pdf",
            out_width = 10, out_height = 4)

## bacterial infection

### setting
forest_plot(forest_calc(dat, "bacterial_infection", "Setting"), xmax = 100,
            out_pdf = "figures/forest_bacterialinfection_setting.pdf",
            out_width = 10, out_height = 5)

### body site
forest_plot(forest_calc(dat, "bacterial_infection", "Body site"), xmax = 100,
            out_pdf = "figures/forest_bacterialinfection_bodysite.pdf",
            out_width = 10, out_height = 5)

## clinical definition
forest_plot(forest_calc(dat, "bacterial_infection", "Clinical definition"), xmax = 100,
            out_pdf = "figures/forest_bacterialinfection_clinicaldef.pdf",
            out_width = 10, out_height = 4)

# sensitivity analysis forest plots

## repeat forest plots for co-infection excluding Type == "Not-specified"
## (in the primary analysis, we assume non-specified infections are co-infections)

### filter dataset
dat_exclude_not_specified <- dat %>%
  filter(!is.na(coinfection) & Type != "Not-specified")

### setting
forest_plot(forest_calc(dat_exclude_not_specified, "coinfection", "Setting"), xmax = 100,
            out_pdf = "figures/forest_coinfection_setting_exclude_not_specified.pdf",
            out_width = 10, out_height = 5)

### body site
forest_plot(forest_calc(dat_exclude_not_specified, "coinfection", "Body site"), xmax = 100,
            out_pdf = "figures/forest_coinfection_bodysite_exclude_not_specified.pdf",
            out_width = 10, out_height = 5)

## quality (risk of bias)
## (for simplicity, just showing this for bacterial infection but could be repeated for other outcomes)

### low
forest_plot(forest_calc(dat %>% filter(bias == "low"),
                        "bacterial_infection", "Setting"), xmax = 100,
            out_pdf = "figures/forest_bacterialinfection_setting_biaslow.pdf",
            out_width = 10, out_height = 5)

### medium
forest_plot(forest_calc(dat %>% filter(bias == "moderate"),
                        "bacterial_infection", "Setting"), xmax = 100,
            out_pdf = "figures/forest_bacterialinfection_setting_biasmoderate.pdf",
            out_width = 10, out_height = 5)

### high
forest_plot(forest_calc(dat %>% filter(bias == "high"),
                        "bacterial_infection", "Setting"), xmax = 100,
            out_pdf = "figures/forest_bacterialinfection_setting_biashigh.pdf",
            out_width = 10, out_height = 5)

## testing (remove unspecified testing and unspecified body site)
## (for simplicity, just showing this for bacterial infection but could be repeated for other outcomes)
forest_plot(forest_calc(dat %>% filter(`Body site` != "Unclear" & Testing != "unspecified"),
                        "bacterial_infection", "Body site"), xmax = 100,
            out_pdf = "figures/forest_bacterialinfection_bodysite_nounspecified.pdf",
            out_width = 10, out_height = 5)

# antibiotic prescribing

## subset studies with antibiotic prescribing information
dat_abx <- dat %>%
  filter(!is.na(abx))

## antibiotic prescribing by end month of study
forest_plot(forest_calc(dat_abx, outcome = "abx", type = "End month", population = "Sample"), xmax = 100,
            out_pdf = "figures/forest_abx.pdf",
            out_width = 10, out_height = 4)
