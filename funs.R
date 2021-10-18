# Functions for: Predictors and Microbiology of Co-Infection in Patients with COVID-19: Living Rapid Review Update and Meta-Regression #
# Langford et al. #
# Script author: Jean-Paul R. Soucy #

## run rma.glmm model for specified outcome w/ or w/out moderators
rma_mod <- function(dat, outcome, term, population = "patients") {
  rma.glmm(
    xi = dat[[outcome]],
    ni = dat[[population]],
    measure = "PLO",
    method = "ML",
    mods = as.formula(paste0("~", paste(term, collapse = " + "))),
    data = dat)
}

## run rma_mod() after filtering dataset to complete cases (no NA)
run_rma_mod <- function(dat, outcomes, var, adjust, population = "patients", verbose = TRUE, adjusted = TRUE) {
  
  ## filter dataset
  if (adjusted) {
    ## filter dataset (use same dataset for adjusted and unadjusted analysis)
    lapply(outcomes, FUN = function(outcome) {
      ## filter to complete cases
      dat_sub <- dat[!is.na(dat[[outcome]]), ] # denominator: those with data for outcome
      n_before <- nrow(dat_sub)
      dat_sub <- dat_sub[, c(outcome, var, adjust, population)] %>%
        {.[complete.cases(.), ]}
      n_after <- nrow(dat_sub)
      
      ## report filtering results
      if (verbose) {
        cat(n_after, "/", n_before, " patients have complete data", sep = "", fill = TRUE)
      }
      
      ## run rma_mod for unadjusted and adjusted analysis
      vars_unadjusted <- c(var)
      vars_adjusted <- unique(c(var, adjust)) # unique in case of overlap between var and adjust
      var_list <- list(vars_unadjusted, vars_adjusted)
      lapply(var_list, function(model_vars) {
        rma_mod(dat = dat_sub, outcome = outcome, term = model_vars, population = population)
      })
    })
  } else {
    ## filter dataset (unadjusted only)
    lapply(outcomes, FUN = function(outcome) {
      ## filter to complete cases
      dat_sub <- dat[!is.na(dat[[outcome]]), ] # denominator: those with data for outcome
      n_before <- nrow(dat_sub)
      dat_sub <- dat_sub[, c(outcome, var, population)] %>%
        {.[complete.cases(.), ]}
      n_after <- nrow(dat_sub)
      
      ## report filtering results
      if (verbose) {
        cat(n_after, "/", n_before, " patients have complete data", sep = "", fill = TRUE)
      }
      
      ## run rma_mod for unadjusted analysis
      rma_mod(dat = dat_sub, outcome = outcome, term = var, population = population)
    })
  }
}

## extract prevalence odds ratio from rma.glmm model
extract_por <- function(var_vals, var) {
  # get var type
  var_type <- get_var_type(var)
  
  # get POR
  if (var_type == "continuous") {
    sprintf("%.2f", exp(var_vals * 10))
  } else if (var_type == "factor") {
    sprintf("%.2f", exp(var_vals))
  } else {
    "VALUE ERROR"
  }
}

## print prevalence with 95% CI
make_prev_ci <- function(mod) {
  paste0(mod$estimate, " (", mod$conf.low, " to ", mod$conf.high, ")")
}

## get variable type
get_var_type <- function(var) {
  case_when(
    var %in% c("setting", "body_site", "study_end_month_grp",
               "region", "risk_of_bias") ~ "factor",
    var %in% c("age", "percent_female", "percent_mechanical_vent",
               "percent_ards", "percent_icu", "percent_severe",
               "percent_smoker", "percent_copd", "percent_cvd",
               "percent_diabetes", "percent_corticosteroid", "percent_il6") ~ "continuous",
    TRUE ~ NA_character_
  )
}

## get outcome display name
get_outcome_name <- function(var) {
  case_when(
    var == "coinfection" ~ "Co-infection",
    var == "secondary_infection" ~ "Secondary infection",
    var == "bacterial_infection" ~ "Bacterial infection",
    TRUE ~ "MISSING OUTCOME NAME"
  )
}

## get variable display name
get_var_name <- function(var) {
  case_when(
    var == "setting" ~ "Setting",
    var == "body_site" ~ "Body site",
    var == "study_end_month_grp" ~ "Study end month",
    var == "region" ~ "Region",
    var == "risk_of_bias" ~ "Risk of bias",
    var == "age" ~ "Age",
    var == "percent_female" ~ "% Female (10% increase)",
    var == "percent_mechanical_vent" ~ "% Mechanical ventilation (10% increase)",
    var == "percent_ards" ~ "% ARDS (10% increase)",
    var == "percent_icu" ~ "% ICU (10% increase)",
    var == "percent_severe" ~ "% Severe (10% increase)",
    var == "percent_smoker" ~ "% Smoker (10% increase)",
    var == "percent_copd" ~ "% COPD (10% increase)",
    var == "percent_cvd" ~ "% CVD (10% increase)",
    var == "percent_diabetes" ~ "% Diabetes (10% increase)",
    var == "percent_corticosteroid" ~ "% Corticosteroid (10% increase)",
    var == "percent_il6" ~ "% IL-6 inhibitor (10% increase)",
    TRUE ~ "MISSING VARIABLE NAME"
  )
}

## get reference level of factor
get_factor_ref_level <- function(var) {
  # grab factor levels from global data frame (dat)
  dat <- get("dat", envir = .GlobalEnv)
  case_when(
    var == "setting" ~ levels(dat[["setting"]])[1],
    var == "body_site" ~ levels(dat[["body_site"]])[1],
    var == "study_end_month_grp" ~ levels(dat[["study_end_month_grp"]])[1],
    var == "region" ~ levels(dat[["region"]])[1],
    var == "risk_of_bias" ~ levels(dat[["risk_of_bias"]])[1],
    TRUE ~ "MISSING REFERENCE LEVEL"
  )
}

## extract data from a list of rma models
extract_rma_mod <- function(mod_list, var) {
  # get var type
  var_type <- get_var_type(var)
  
  # extract data
  lapply(mod_list, function(mod) {
    tidy(mod, conf.int = TRUE, measure = "PLO", exponentiate = FALSE) %>%
      select(term, estimate, conf.low, conf.high) %>%
      filter(grepl(paste0("^", var), term)) %>%
      mutate(
        term = sub(paste0("^", var), "", term),
        estimate = extract_por(estimate, var),
        conf.low = extract_por(conf.low, var),
        conf.high = extract_por(conf.high, var)
        )
  })
}

## make table row
make_table_row <- function(results_list, var) {
  
  # get variable type
  var_type <- get_var_type(var)
  match.arg(var_type,
            choices = c("factor", "continuous"),
            several.ok = FALSE)
  
  # make table row
  if (var_type == "factor") {
    term_col <- c(get_var_name(var), get_factor_ref_level(var), results_list[[1]]$term)
    value_cols <- lapply(results_list, function(result) {
      c("", "Reference", make_prev_ci(result))
    })
  } else {
    term_col <- get_var_name(var)
    value_cols <- lapply(results_list, function(result) {
      make_prev_ci(result)
    })
  }
  data.frame(
    terms = term_col,
    matrix(unlist(value_cols), ncol = length(value_cols), byrow = FALSE))
}

## make sample size columns
make_sample_columns <- function(dat, outcomes, var, adjust, population = "patients", adjusted = TRUE) {
  # grab sample sizes from global data frame (dat)
  dat_sample <- get("dat", envir = .GlobalEnv)
  # add sample size columns
  if (adjusted) {
    for (i in 1:length(outcomes)) {
      dat_i <- dat_sample[, c(outcomes[i], var, adjust, population)] %>%
        {.[complete.cases(.), ]}
      sample_size <- nrow(dat_i)
      if (get_var_type(var) == "factor") {
        factor_levels_n <- as.vector(table(dat_i[, var]))
        sample_size_col <- c(sample_size, factor_levels_n)
      } else {
        sample_size_col <- sample_size
      }
      dat <- dat %>%
        add_column(sample_size_col, .after = 3 * i, .name_repair = "unique")
    }
  } else {
    for (i in 1:length(outcomes)) {
      dat_i <- dat_sample[, c(outcomes[i], var, population)] %>%
        {.[complete.cases(.), ]}
      sample_size <- nrow(dat_i)
      if (get_var_type(var) == "factor") {
        factor_levels_n <- as.vector(table(dat_i[, var]))
        sample_size_col <- c(sample_size, factor_levels_n)
      } else {
        sample_size_col <- sample_size
      }
      dat <- dat %>%
        add_column(sample_size_col, .after = 2 * i, .name_repair = "unique")
    }
  }
  
  # return data
  return(dat)
}

## make table header
make_table_headers <- function(tab, outcomes, adjusted = TRUE, exclude_not_specified = FALSE) {
  # grab sample sizes from global data frame (dat)
  dat_sample <- get("dat", envir = .GlobalEnv)
  if (exclude_not_specified) {
    dat_sample <- dat_sample %>% filter(Type != "Not-specified")
  }
  # create header rows
  if (adjusted) {
    # create row 1
    row_1 <- c("Characteristic", rep("", each = 3 * length(outcomes)))
    # calculate sample size for each outcome
    for (i in 1:length(outcomes)) {
      outcome_n <- nrow(dat_sample[!is.na(dat_sample[[outcomes[i]]]), ])
      outcome_header <- paste0(get_outcome_name(outcomes[i]), " (n = ", outcome_n, ")")
      row_1[2 + 3 * (i - 1)] <- outcome_header
    }
    row_2 <- c("", c(rep(c("Unadjusted", "Adjusted", "Studies included"), times = length(outcomes))))
    header <- matrix(c(row_1, row_2), nrow = 2, byrow = TRUE,
                     dimnames = list(1:2, names(tab)))
  } else {
    # create row 1
    row_1 <- c("Characteristic", rep("", each = 2 * length(outcomes)))
    # calculate sample size for each outcome
    for (i in 1:length(outcomes)) {
      outcome_n <- nrow(dat_sample[!is.na(dat_sample[[outcomes[i]]]), ])
      outcome_header <- paste0(get_outcome_name(outcomes[i]), " (n = ", outcome_n, ")")
      row_1[2 + 2 * (i - 1)] <- outcome_header
    }
    row_2 <- c("", c(rep(c("Unadjusted", "Studies included"), times = length(outcomes))))
    header <- matrix(c(row_1, row_2), nrow = 2, byrow = TRUE,
                     dimnames = list(1:2, names(tab)))
  }
  
  ## return data with headers
  return(rbind(header, tab))
}

## replace rows in the summary table
replace_rows <- function(dat, outcomes, var, adjust, summary_table, characteristic, n_rows) {
  # generate new rows
  rows <- run_rma_mod(dat, outcomes, var, adjust) %>%
    unlist(recursive = FALSE) %>%
    extract_rma_mod(var) %>%
    make_table_row(var) %>%
    make_sample_columns(outcomes, var, adjust)
  
  # find row number to begin replacement
  n_rows_begin <- grep(characteristic, summary_table$terms)
  
  # replace rows
  summary_table[n_rows_begin:(n_rows_begin + n_rows - 1), ] <- rows
  
  # return summary table
  return(summary_table)
}

## calculate forest plots for specified outcome w/ or w/out subgroups
forest_calc <- function(dat, outcome, type, population = "patients") {
  
  if (type == "All") {
    metaprop(
      event = dat[[outcome]],
      n = dat[[population]],
      studlab = dat[["study"]],
      method = "GLMM",
      sm = "PLOGIT"
    )
  } else {
    metaprop(
      event = dat[[outcome]],
      n = dat[[population]],
      studlab = dat[["study"]],
      method = "GLMM",
      sm = "PLOGIT",
      byvar = dat[[type]],
      bylab = type
    )
  }
}

## plot forest plot
forest_plot <- function(dat,
                        xmin = 0, # x-axis minimum value
                        xmax = 100, # x-axis maximum value
                        order_subgroups = FALSE, # order subgroups from low to high prevalence
                        show_ind_studies = FALSE, # show individual studies?
                        out_pdf = NULL, # if specified, path to output pdf
                        out_width = 7, # width of output
                        out_height = 7 # height of output
                        ){
  
  ## determine subgroup variable (if null, show all studies)
  type <- ifelse(is.null(dat$byvar), "All", dat$byvar)
  
  ## order subgroup results in decreasing order of prevalence
  if (order_subgroups == TRUE) {
    o <- order(dat$TE.random.w, decreasing = FALSE)
    for (var in c("bylevs", grep("\\.w$", names(dat)[!names(dat) %in% c("df.hakn.w", "df.Q.w")], value = TRUE))) {
      dat[[var]] <- dat[[var]][o]
    }
  }
  
  ## open graphics device to save plot
  if (!is.null(out_pdf)) {
    pdf(file = out_pdf, width = out_width, height = out_height)
  }
  
  ## forest plot
  meta::forest(dat,
               xlim = c(xmin, xmax),
               pscale = 100,
               rightcols = FALSE,
               leftcols = c("studlab", "n", "effect", "ci"),
               leftlabs = c("Subgroup",
                            "Total Patients",
                            "Prevalence (%)",
                            "95% C.I."),
               xlab = "Prevalence (%)", smlab = "",
               weight.study = "random", squaresize = 0.5, col.square = "navy",
               col.square.lines = "navy",
               col.diamond = "maroon", col.diamond.lines = "maroon",
               pooled.totals = TRUE,
               comb.fixed = FALSE,
               fs.hetstat = 10,
               print.tau2 = TRUE,
               print.Q = TRUE,
               print.pval.Q = TRUE,
               print.I2 = TRUE,
               digits = 1,
               # bylab = "",
               bylab = "",
               col.by = "black",
               study.results = show_ind_studies
  )
  
  ## close graphics device
  if (!is.null(out_pdf)) {
    dev.off()
  }
  
}
