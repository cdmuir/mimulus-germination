##### Write Stan blocks for germination model ##################################

# Function to create custom generated quantities blocks in Stan for each subcohort-pop combination
make_germ_generated_quantities <- function(pop, sub_cohort, start, model, 
                                           nonadd) {
  
  # Checks
  check_pop(pop)
  model %<>% check_model()
  
  # Inverse link function
  inverse_link <- switch(model, 
    discrete_loglogistic = "exp", 
    discrete_lognormal = ""
  )
  
  # Suffix for labeling various things in Stan
  suffix <- make_suffix(pop, sub_cohort)
  
  # Additional model parameters
  parameters <- add_parameters(model, pop)
  
  generated_quantities_block <- sprintf("
                         
  for (i in 1:n%s) {
    mu%s_germ[i] = %s(bPop_germ[%i] + bCohortNorth_germ * Cohort%s[i] +
      bBlock_germ[Block%s[i]] +
      bGeno%s_germ[Sire%s[i]] + bGeno%s_germ[Dam%s[i]]);", 
                                        suffix, 
                                        suffix, inverse_link, which(pop == pop_levels()), suffix, 
                                        suffix,
                                        pop, suffix, pop, suffix
  )
  
  if (nonadd) {
    
    generated_quantities_block %<>% 
      stringr::str_replace("\\);$", " +\n") %>%
      stringr::str_c(sprintf("      bSireDam%s_germ[SireDam%s[i]]);", pop, suffix))
    
  }
  
  generated_quantities_block %<>% stringr::str_c(sprintf("
    log_lik_germ[%si] = %s_%s_lpmf(DaysToGerm%s[i] | mu%s_germ[i]%s);
    predict_germ[%si] = %s_%s_rng(mu%s_germ[i]%s);
  }", start, model, sub_cohort, suffix, suffix, parameters,
                                                         start, model, sub_cohort, suffix, parameters))
  
  generated_quantities_block
  
}

# Function to create custom model blocks in Stan for each subcohort-pop combination
make_germ_model <- function(pop, sub_cohort, model, nonadd) {
  
  # Checks
  check_pop(pop)
  model %<>% check_model()
  
  # Inverse link function
  inverse_link <- switch(model, 
                         discrete_loglogistic = "exp", 
                         discrete_lognormal = ""
  )
  
  # Suffix for labeling various things in Stan
  suffix <- make_suffix(pop, sub_cohort)
  
  # Additional model parameters
  parameters <- add_parameters(model, pop)
  
  model_block <- sprintf("
  for (i in 1:n%s) {
    mu%s[i] = %s(bPop_germ[%i] + bCohortNorth_germ * Cohort%s[i] +
      bBlock_germ[Block%s[i]] +
      bGeno%s_germ[Sire%s[i]] + bGeno%s_germ[Dam%s[i]]);", 
                         suffix, 
                         suffix, inverse_link, which(pop == pop_levels()), suffix, 
                         suffix, 
                         pop, suffix, pop, suffix
  )
  if (nonadd) {
    
    model_block %<>% 
      stringr::str_replace("\\);", " +\n") %>%
      stringr::str_c(sprintf("      bSireDam%s_germ[SireDam%s[i]]);", pop, suffix))
    
  }
  
  model_block %<>% 
    stringr::str_c(
      sprintf("\n\t\ttarget += %s_%s_lpmf(DaysToGerm%s[i] | mu%s[i]%s);
  }\n", model, sub_cohort, suffix, suffix, parameters)
    )
  
  model_block
  
}

# Function to create custom likelihood functions in Stan for each subcohort
make_germ_function <- function(model, date_sown, dates_censused, 
                               sub_cohort = "") {
  
  # Check model
  model %<>% check_model()
  
  # Additional model parameters
  parameters <- add_parameters(model, pop = NA)
  types <- add_parameter_types(model)
  
  # Convert to date class and sort
  date_sown %<>% as.Date()
  dates_censused %<>% as.Date()
  dates_censused %<>% subtract(date_sown) %>% 
    as.integer() %>% 
    sort()
  n <- length(dates_censused)
  
  # Initiate blocks
  # lpmf_function calculates log-posterior mass function
  # rng_function generates random numbers for prediction
  lpmf_function <- sprintf("
  real %s_%s_lpmf(int y, real mu%s) {
    
    real ret;

    if (y == %i) {
      ret = %s_lcdf(%i | mu%s);
    }
  ", model, sub_cohort, types,
                           dates_censused[1], 
                           model, dates_censused[1], parameters)
  
  rng_function <- sprintf("
  real %s_%s_rng(real mu%s) {
                           
    real ret;
    ret = %s_rng(mu%s);
                           
    if (ret <= %i) {
      ret = %i;
    }
  ", model, sub_cohort, types,
                          model, parameters,
                          dates_censused[1], 
                          dates_censused[1])
  
  for (i in 2L:n) {
    
    if (dates_censused[i] - dates_censused[i - 1L] == 1L) {
      
      lpmf_function %<>% stringr::str_c(sprintf("
    else if (y == %i) {
      ret = %s_lpmf(%i | mu%s);
    }", dates_censused[i], 
                                                model, dates_censused[i], parameters))
      
      rng_function %<>% stringr::str_c(sprintf("
    else if (ret == %i) {
      ret = %i;
    }", dates_censused[i], 
                                               dates_censused[i]))
      
    } else {
      
      lpmf_function %<>% stringr::str_c(sprintf("
    else if (y == %i) {
      ret = log_diff_exp(%s_lcdf(%i | mu%s),
                         %s_lcdf(%i | mu%s));
    }", dates_censused[i], 
                                                model, dates_censused[i], parameters,
                                                model, dates_censused[i - 1L], parameters))
      
      rng_function %<>% stringr::str_c(sprintf("
    else if (ret <= %i) {
      ret = %i;
    }", dates_censused[i], 
                                               dates_censused[i]))
      
    } 
  }
  
  lpmf_function %<>% stringr::str_c(sprintf("
    else {
      ret = %s_lccdf(%i | mu%s);
    }

  return ret;

  }
", model, dates_censused[n], parameters))
  
  # if predicted DaysToGerm is larger than maximum observed, censor to maximum value
  rng_function %<>% stringr::str_c(sprintf("
    else {
      ret = %i;
    }

  return ret;

  }
", dates_censused[n]))
  
  # return
  stringr::str_c(lpmf_function, rng_function)
  
}

# Function to create custom data blocks in Stan for each subcohort-pop combination
make_germ_data <- function(pop, sub_cohort, nonadd) {
  
  # Checks
  check_pop(pop)
  
  # Suffix for labeling various things in Stan
  suffix <- make_suffix(pop, sub_cohort)
  
  data_block <- sprintf("

  int<lower=0> n%s; // Length of data in SubCohort-pop combination
  int<lower=0,upper=1> Cohort%s[n%s]; // Logical. In North Cohort?
  int<lower=0> Block%s[n%s];
  int<lower=0> Sire%s[n%s];
  int<lower=0> Dam%s[n%s];
  int<lower=0> DaysToGerm%s[n%s];", 
                        suffix, suffix, suffix, suffix, suffix, suffix, suffix, suffix, suffix, suffix, suffix)
  
  if (nonadd) {
    data_block %<>% stringr::str_c(sprintf("
  int<lower=0> SireDam%s[n%s];", suffix, suffix))
  }
  
  data_block
  
}

##### Write Stan blocks for survival model #####################################

# Function to create custom generated quantities blocks in Stan for each population
make_surv_generated_quantities <- function(pop, start, nonadd, interaction) {
  
  # Checks
  check_pop(pop)
  
  generated_quantities_block <- sprintf("
                                        
  for (i in 1:n%s) {
    theta%s[i] = inv_logit(bPop_surv[%i] + bCohortNorth_surv * Cohort%s[i] +
                           bBlock_surv[Block%s[i]] +
                           bGeno%s_surv[Sire%s[i]] + 
                           bGeno%s_surv[Dam%s[i]]);", 
                                        pop, 
                                        pop, which(pop == pop_levels()), pop, 
                                        pop,
                                        pop, pop, pop, pop)
  
  if (interaction) {
    
    generated_quantities_block %<>%
      stringr::str_replace(
        "bCohortNorth_surv", stringr::str_c("bCohortNorth_surv[", 
                                       which(pop == pop_levels()), "]")
      )
    
    generated_quantities_block %<>% 
      stringr::str_replace_all(sprintf("([\\s]*)(bGeno%s_surv\\[)([a-zA-Z]+%s\\[i\\]\\])([\\s]*)", pop, pop),
                               sprintf("\\1\\2\\3[1] * abs(Cohort%s[i] - 1) + \\2\\3[2] * Cohort%s[i]\\4",
                                       pop, pop))
    
  }
  
  if (nonadd) {
    
    generated_quantities_block %<>% 
      stringr::str_replace("\\);$", " +\n") %>%
      stringr::str_c(sprintf("                           bSireDam%s_surv[SireDam%s[i]]);", pop, pop))
    
    if (interaction) {
      
      generated_quantities_block %<>% 
        stringr::str_replace_all(sprintf("([\\s]*)(bSireDam%s\\[)(SireDam%s\\[i\\]\\])([\\s]*)", pop, pop),
                                 sprintf("\\1\\2\\3[1] * abs(Cohort%s[i] - 1) + \\2\\3[2] * Cohort%s[i]\\4",
                                         pop, pop))
      
    }
    
  }
  
  generated_quantities_block %<>% stringr::str_c(sprintf("
  log_lik_surv[%si] = binomial_lpmf(WinterSurv%s[i] | 1, theta%s[i]);
  predict_surv[%si] = binomial_rng(1, theta%s[i]);
}", start, pop, pop,
                                                         start, pop))
  
  generated_quantities_block
  
}

# Function to create custom model blocks in Stan for each population
make_surv_model <- function(pop, nonadd, interaction) {
  
  # Checks
  check_pop(pop)
  
  model_block <- sprintf("
  for (i in 1:n%s) {
    theta%s[i] = inv_logit(bPop_surv[%i] + bCohortNorth_surv * Cohort%s[i] +
                           bBlock_surv[Block%s[i]] +
                           bGeno%s_surv[Sire%s[i]] + 
                           bGeno%s_surv[Dam%s[i]]);", 
                         pop, 
                         pop, which(pop == pop_levels()), pop, 
                         pop, 
                         pop, pop, pop, pop)
  
  if (interaction) {
    
    model_block %<>%
      stringr::str_replace("bCohortNorth_surv", stringr::str_c("bCohortNorth_surv[", which(pop == pop_levels()), "]"))
    
    model_block %<>% 
      stringr::str_replace_all(sprintf("([\\s]*)(bGeno%s_surv\\[)([a-zA-Z]+%s\\[i\\]\\])([\\s]*)", pop, pop),
                               sprintf("\\1\\2\\3[1] * abs(Cohort%s[i] - 1) + \\2\\3[2] * Cohort%s[i]\\4",
                                       pop, pop))
    
  }
  
  if (nonadd) {
    
    model_block %<>% 
      stringr::str_replace("\\);", " +\n") %>%
      stringr::str_c(sprintf("                           bSireDam%s_surv[SireDam%s[i]]);", pop, pop))
    
    if (interaction) {
      
      model_block %<>% 
        stringr::str_replace_all(sprintf("([\\s]*)(bSireDam%s_surv\\[)(SireDam%s\\[i\\]\\])([\\s]*)", pop, pop),
                                 sprintf("\\1\\2\\3[1] * abs(Cohort%s[i] - 1) + \\2\\3[2] * Cohort%s[i]\\4",
                                         pop, pop))
      
    }
  }
  
  model_block %<>% 
    stringr::str_c(sprintf("\n\t\ttarget += binomial_lpmf(WinterSurv%s[i] | 1, theta%s[i]);\n\t}\n", 
                           pop, pop))
  
  model_block
  
}

# Function to create custom data blocks in Stan for each subcohort-pop combination
make_surv_data <- function(pop, nonadd) {
  
  # Checks
  check_pop(pop)
  
  data_block <- sprintf("
                        
  int<lower=0> n%s; // Length of data in population
  int<lower=0,upper=1> Cohort%s[n%s]; // Logical. In North Cohort?
  int<lower=0> Block%s[n%s];
  int<lower=0> Sire%s[n%s];
  int<lower=0> Dam%s[n%s];
  int<lower=0> WinterSurv%s[n%s];", 
                        pop, pop, pop, pop, pop, pop, pop, pop, pop, pop, pop)
  
  if (nonadd) {
    data_block %<>% stringr::str_c(sprintf("
  int<lower=0> SireDam%s[n%s];", pop, pop))
  }
  
  data_block
  
}


##### Combined model ##########################################################

# Write Stan model ----
write_model <- function(
  model, sow_dates, census_dates, data, file, nonadd = TRUE,
  poph2 = TRUE, interaction = TRUE
) {
  
  ## Functions Block (germination only) -----
  
  functions_block <- stringr::str_c(
    "functions {\n",
    readChar("stan/discrete_lognormal.stan", 1e5), 
    readChar("stan/loglogistic.stan", 1e5),
    readChar("stan/discrete_loglogistic.stan", 1e5), "\n")

  for (i in 1:nrow(sow_dates)) {
    
    sub_cohort <- sow_dates$sub_cohort[i]
    sow_date <- sow_dates$date[i]
    cohort <- ifelse(stringr::str_detect(sub_cohort, "north"), "north", "south")
    
    functions_block %<>% 
      stringr::str_c(
        make_germ_function(
          model = model, 
          date_sown = sow_date, 
          dates_censused = 
            dplyr::filter(census_dates, Cohort == cohort)$GermDate,
          sub_cohort = sub_cohort
        )
      )
    
  }
  
  functions_block %<>% stringr::str_c("\n}")
  
  ## Data Block (germination) -----
  
  data_block_germ <- " data {\n\n// Germination model\n\n\t// Data lengths\n"
  
  for (j in levels(data$Pop)) {
    data_block_germ %<>% stringr::str_c("\tint<lower=1> nSire", j, "_germ;\n")
    if (nonadd) {
      data_block_germ %<>% stringr::str_c("\tint<lower=1> nSireDam", j, "_germ;\n")
    }
  }
  
  data_block_germ %<>% 
    stringr::str_c("\tint<lower=1> nPop_germ;\n") %>%
    stringr::str_c("\tint<lower=1> nBlock_germ;") %>%
    stringr::str_c("\treal ntd[nPop_germ];")
  
  
  for (i in sow_dates$sub_cohort) {
    for (j in levels(data$Pop)) {
      data_block_germ %<>% stringr::str_c(
        make_germ_data(pop = j, sub_cohort = i, nonadd = nonadd)
      )
    }
  }

  ## Data Block (survival) -----
  
  data_block_surv <- "\t// Survival model\n\t// Data lengths\n"
  
  for (pop in levels(data$Pop)) {
    data_block_surv %<>% stringr::str_c("\tint<lower=1> nSire", pop, "_surv;\n")
    if (nonadd) {
      data_block_surv %<>% stringr::str_c("\tint<lower=1> nSireDam", pop, "_surv;\n")
    }
  }
  
  data_block_surv %<>% 
    stringr::str_c("\tint<lower=1> nPop_surv;\n") %>%
    stringr::str_c("\tint<lower=1> nBlock_surv;")
  
  for (pop in levels(data$Pop)) {
    data_block_surv %<>% stringr::str_c(make_surv_data(pop = pop, nonadd = nonadd))
  }
  
  data_block <- stringr::str_c(data_block_germ, "\n\n", data_block_surv, "\n\n}")
  
  ## Parameters Block (germination) -----
  
  parameters_block_germ <- " parameters {
  
  // Germination model
  
  // Fixed effect of cohort
    real bCohortNorth_germ;
  
  // Fixed effect of population
    vector[nPop_germ] bPop_germ;
  
  // Random effect of block
    vector[nBlock_germ] bBlock_germ;
    real<lower=0> sBlock_germ;
  
  // Additive genetic variance
    vector[nSireCUR_germ] bGenoCUR_germ;
    vector[nSireWFM_germ] bGenoWFM_germ;
    vector[nSireNMT_germ] bGenoNMT_germ;
    vector[nSireLIJ_germ] bGenoLIJ_germ;
    vector[nSireRCK_germ] bGenoRCK_germ;
  
    real<lower=0> sGeno_germ;\n"
  
  if (poph2 & !nonadd) {
    
    parameters_block_germ %<>% 
      stringr::str_replace("real<lower=0> sGeno_germ;",
                           "real<lower=0> sGenoCUR_germ;
    real<lower=0> sGenoWFM_germ;
    real<lower=0> sGenoNMT_germ;
    real<lower=0> sGenoLIJ_germ;
    real<lower=0> sGenoRCK_germ;\n") 
    
  }
  
  if (!poph2 & nonadd) {
    
    parameters_block_germ %<>% 
      stringr::str_c("

  // Nonadditive genetic variance
    vector[nSireDamCUR_germ] bSireDamCUR_germ;
    vector[nSireDamWFM_germ] bSireDamWFM_germ;
    vector[nSireDamNMT_germ] bSireDamNMT_germ;
    vector[nSireDamLIJ_germ] bSireDamLIJ_germ;
    vector[nSireDamRCK_germ] bSireDamRCK_germ;
    
    real<lower=0> sSireDam_germ;\n")
    
  }
  
  if (poph2 & nonadd) {
    
    parameters_block_germ %<>% 
      stringr::str_replace("real<lower=0> sGeno_germ;",
                           "real<lower=0> sGenoCUR_germ;
    real<lower=0> sGenoWFM_germ;
    real<lower=0> sGenoNMT_germ;
    real<lower=0> sGenoLIJ_germ;
    real<lower=0> sGenoRCK_germ;") %>% 
      
      stringr::str_c("

  // Nonadditive genetic variance
    vector[nSireDamCUR_germ] bSireDamCUR_germ;
    vector[nSireDamWFM_germ] bSireDamWFM_germ;
    vector[nSireDamNMT_germ] bSireDamNMT_germ;
    vector[nSireDamLIJ_germ] bSireDamLIJ_germ;
    vector[nSireDamRCK_germ] bSireDamRCK_germ;
    
    real<lower=0> sSireDamCUR_germ;
    real<lower=0> sSireDamWFM_germ;
    real<lower=0> sSireDamNMT_germ;
    real<lower=0> sSireDamLIJ_germ;
    real<lower=0> sSireDamRCK_germ;\n")
    
  }
  
  parameters_block_germ %<>% 
    stringr::str_c(add_parameter_definitions(model))
  
  ## Parameters Block (survival) ----

  parameters_block_surv <- "\t// Survival model
  
  // Fixed effect of cohort
    real bCohortNorth_surv;

  // Fixed effect of population
    vector[nPop_surv] bPop_surv;

  // Random effect of block
    vector[nBlock_surv] bBlock_surv;
    real<lower=0> sBlock_surv;
  
  // Additive genetic variance
    vector[nSireCUR_surv] bGenoCUR_surv;
    vector[nSireWFM_surv] bGenoWFM_surv;
    vector[nSireNMT_surv] bGenoNMT_surv;
    vector[nSireLIJ_surv] bGenoLIJ_surv;
    vector[nSireRCK_surv] bGenoRCK_surv;
  
    real<lower=0> sGeno_surv;\n"

  if (poph2 & !nonadd) {
    
    parameters_block_surv %<>% 
      stringr::str_replace("real<lower=0> sGeno_surv;",
                           "real<lower=0> sGenoCUR_surv;
    real<lower=0> sGenoWFM_surv;
    real<lower=0> sGenoNMT_surv;
    real<lower=0> sGenoLIJ_surv;
    real<lower=0> sGenoRCK_surv;\n") 
    
  }
  
  if (!poph2 & nonadd) {
    
    parameters_block_surv %<>% 
      stringr::str_c("

  // Nonadditive genetic variance
    vector[nSireDamCUR_surv] bSireDamCUR_surv;
    vector[nSireDamWFM_surv] bSireDamWFM_surv;
    vector[nSireDamNMT_surv] bSireDamNMT_surv;
    vector[nSireDamLIJ_surv] bSireDamLIJ_surv;
    vector[nSireDamRCK_surv] bSireDamRCK_surv;
    
    real<lower=0> sSireDam_surv;\n")
    
  }
  
  if (poph2 & nonadd) {
    
    parameters_block_surv %<>% 
      stringr::str_replace("real<lower=0> sGeno_surv;",
                           "real<lower=0> sGenoCUR_surv;
    real<lower=0> sGenoWFM_surv;
    real<lower=0> sGenoNMT_surv;
    real<lower=0> sGenoLIJ_surv;
    real<lower=0> sGenoRCK_surv;") %>% 
      
      stringr::str_c("

  // Nonadditive genetic variance
    vector[nSireDamCUR_surv] bSireDamCUR_surv;
    vector[nSireDamWFM_surv] bSireDamWFM_surv;
    vector[nSireDamNMT_surv] bSireDamNMT_surv;
    vector[nSireDamLIJ_surv] bSireDamLIJ_surv;
    vector[nSireDamRCK_surv] bSireDamRCK_surv;
    
    real<lower=0> sSireDamCUR_surv;
    real<lower=0> sSireDamWFM_surv;
    real<lower=0> sSireDamNMT_surv;
    real<lower=0> sSireDamLIJ_surv;
    real<lower=0> sSireDamRCK_surv;\n")
    
  }
  
  # Add parameters to fit genotype x environment interaction
  if (interaction) {
    
    # Population effects vary with garden location
    parameters_block_surv %<>% 
      stringr::str_replace_all("real bCohortNorth_surv;",
                               "vector[nPop_surv] bCohortNorth_surv;") %>%
      
      # Add south- and north-specific coefficients
      stringr::str_replace_all(
        "vector\\[nSire([A-Z]{3})_surv\\] bGeno([A-Z]{3})_surv;",
        "vector\\[2] bGeno\\2_surv\\[nSire\\1_surv\\];"
      ) %>%
      stringr::str_replace_all(
        "vector\\[nSireDam([A-Z]{3})_surv\\] bSireDam([A-Z]{3})_surv;",
        "vector\\[2] bSireDam\\2_surv\\[nSireDam\\1_surv\\];"
      ) %>%
      
      # Add tau - vector of separate variances for south and north
      stringr::str_replace_all(
        "real<lower=0> sGeno([A-Z]{0,3})_surv;", 
        "vector<lower=0>[2] tauGeno\\1_surv;"
      ) %>%
      stringr::str_replace_all(
        "real<lower=0> sSireDam([A-Z]{0,3})_surv;", 
        "vector<lower=0>[2] tauSireDam\\1_surv;"
      ) %>%
      
      # Add covariance matrix between south and north coefficients
      stringr::str_c(stringr::str_c(
        glue::glue("\n\t\tcorr_matrix[2] Omega{x};",
                   x = tidyr::crossing(
                     x = if (nonadd) {c("Geno", "SireDam")} else {"Geno"},
                     y = if (poph2) {pop_levels()} else {""}
                   ) %>%
                     dplyr::transmute(z = paste0(x, y)) %>%
                     dplyr::pull(z)
        ), 
        collapse = "\n"))
    
  }
  
  parameters_block <- stringr::str_c(parameters_block_germ, "\n\n",
                                     parameters_block_surv, "\n}")
  
  ## Model Block (germination priors) -----
  
  model_priors_germ <- " model {\n\n\t// Germination model priors\n"
  
  for (i in sow_dates$sub_cohort) {
    for (j in levels(data$Pop)) {
      suffix <- stringr::str_c("_", i, "_", j)
      model_priors_germ %<>% 
        stringr::str_c("\tvector[n", suffix, "] mu",suffix, ";\n")
    }
  }
  model_priors_germ %<>% stringr::str_c("
  // Priors on block
  sBlock_germ ~ cauchy(0, 0.1);
  bBlock_germ ~ normal(0, sBlock_germ);
    
  // Priors on population
  bPop_germ ~ normal(0, 2);
    
  // Prior on cohort
  bCohortNorth_germ ~ normal(0, 1);
                                        
  // Priors on genetic variance\n")
  
  if (!poph2 & !nonadd) {
    
    model_priors_germ %<>%
      stringr::str_c("
  sGeno_germ ~ cauchy(0, 0.1);
    
  bGenoCUR_germ ~ normal(0, sGeno_germ);
  bGenoWFM_germ ~ normal(0, sGeno_germ);
  bGenoNMT_germ ~ normal(0, sGeno_germ);
  bGenoLIJ_germ ~ normal(0, sGeno_germ);
  bGenoRCK_germ ~ normal(0, sGeno_germ);
    ")
    
  }
  
  if (poph2 & !nonadd) {
    
    model_priors_germ %<>%
      stringr::str_c("
  sGenoCUR_germ ~ cauchy(0, 0.1);
  sGenoWFM_germ ~ cauchy(0, 0.1);
  sGenoNMT_germ ~ cauchy(0, 0.1);
  sGenoLIJ_germ ~ cauchy(0, 0.1);
  sGenoRCK_germ ~ cauchy(0, 0.1);
  
  bGenoCUR_germ ~ normal(0, sGenoCUR_germ);
  bGenoWFM_germ ~ normal(0, sGenoWFM_germ);
  bGenoNMT_germ ~ normal(0, sGenoNMT_germ);
  bGenoLIJ_germ ~ normal(0, sGenoLIJ_germ);
  bGenoRCK_germ ~ normal(0, sGenoRCK_germ);
    ")
    
  }
  
  if (!poph2 & nonadd) {
    
    model_priors_germ %<>% 
      stringr::str_c("

  sGeno_germ ~ cauchy(0, 0.1);
  sSireDam_germ ~ cauchy(0, 0.1);
  
  bGenoCUR_germ ~ normal(0, sGeno_germ);
  bGenoWFM_germ ~ normal(0, sGeno_germ);
  bGenoNMT_germ ~ normal(0, sGeno_germ);
  bGenoLIJ_germ ~ normal(0, sGeno_germ);
  bGenoRCK_germ ~ normal(0, sGeno_germ);
  
  bSireDamCUR_germ ~ normal(0, sSireDam_germ);
  bSireDamWFM_germ ~ normal(0, sSireDam_germ);
  bSireDamNMT_germ ~ normal(0, sSireDam_germ);
  bSireDamLIJ_germ ~ normal(0, sSireDam_germ);
  bSireDamRCK_germ ~ normal(0, sSireDam_germ);")
    
  }
  
  if (poph2 & nonadd) {
    
    model_priors_germ %<>% 
      stringr::str_c("

  sGenoCUR_germ ~ cauchy(0, 0.1);
  sGenoWFM_germ ~ cauchy(0, 0.1);
  sGenoNMT_germ ~ cauchy(0, 0.1);
  sGenoLIJ_germ ~ cauchy(0, 0.1);
  sGenoRCK_germ ~ cauchy(0, 0.1);
  
  sSireDamCUR_germ ~ cauchy(0, 0.1);
  sSireDamWFM_germ ~ cauchy(0, 0.1);
  sSireDamNMT_germ ~ cauchy(0, 0.1);
  sSireDamLIJ_germ ~ cauchy(0, 0.1);
  sSireDamRCK_germ ~ cauchy(0, 0.1);
  
  bGenoCUR_germ ~ normal(0, sGenoCUR_germ);
  bGenoWFM_germ ~ normal(0, sGenoWFM_germ);
  bGenoNMT_germ ~ normal(0, sGenoNMT_germ);
  bGenoLIJ_germ ~ normal(0, sGenoLIJ_germ);
  bGenoRCK_germ ~ normal(0, sGenoRCK_germ);
  
  bSireDamCUR_germ ~ normal(0, sSireDamCUR_germ);
  bSireDamWFM_germ ~ normal(0, sSireDamWFM_germ);
  bSireDamNMT_germ ~ normal(0, sSireDamNMT_germ);
  bSireDamLIJ_germ ~ normal(0, sSireDamLIJ_germ);
  bSireDamRCK_germ ~ normal(0, sSireDamRCK_germ);")
    
  }
  
  model_priors_germ %<>% stringr::str_c("\n", add_parameter_priors(model))
  
  ## Model Block (survival priors) ----

  model_priors_surv <- "\n\t// Survival model priors\n"
  
  for (pop in levels(data$Pop)) {
    model_priors_surv %<>% stringr::str_c("\tvector[n", pop, "] theta", pop, ";\n")
  }
  
  if (interaction) model_priors_surv %<>% 
    stringr::str_c("\tvector[2] Zero;\n\tZero[1] = 0;\n\tZero[2] = 0;\n")
  
  model_priors_surv %<>% stringr::str_c("
  // Priors on block
  sBlock_surv ~ cauchy(0, 0.1);
  bBlock_surv ~ normal(0, sBlock_surv);
  
  // Priors on population
  bPop_surv ~ normal(0, 2);
  
  // Prior on cohort
  bCohortNorth_surv ~ normal(0, 1);
                                        
  // Priors on additive genetic variance")
  
  if (!poph2 & !nonadd) {
    
    model_priors_surv %<>%
      stringr::str_c("
  sGeno_surv ~ cauchy(0, 0.1);
    
  bGenoCUR_surv ~ normal(0, sGeno_surv);
  bGenoWFM_surv ~ normal(0, sGeno_surv);
  bGenoNMT_surv ~ normal(0, sGeno_surv);
  bGenoLIJ_surv ~ normal(0, sGeno_surv);
  bGenoRCK_surv ~ normal(0, sGeno_surv);
    ")
    
  }
  
  if (poph2 & !nonadd) {
    
    model_priors_surv %<>%
      stringr::str_c("
  sGenoCUR_surv ~ cauchy(0, 0.1);
  sGenoWFM_surv ~ cauchy(0, 0.1);
  sGenoNMT_surv ~ cauchy(0, 0.1);
  sGenoLIJ_surv ~ cauchy(0, 0.1);
  sGenoRCK_surv ~ cauchy(0, 0.1);
  
  bGenoCUR_surv ~ normal(0, sGenoCUR_surv);
  bGenoWFM_surv ~ normal(0, sGenoWFM_surv);
  bGenoNMT_surv ~ normal(0, sGenoNMT_surv);
  bGenoLIJ_surv ~ normal(0, sGenoLIJ_surv);
  bGenoRCK_surv ~ normal(0, sGenoRCK_surv);
    ")
    
  }
  
  if (!poph2 & nonadd) {
    
    model_priors_surv %<>% 
      stringr::str_c("

  sGeno_surv ~ cauchy(0, 0.1);
  sSireDam_surv ~ cauchy(0, 0.1);

  bGenoCUR_surv ~ normal(0, sGeno_surv);
  bGenoWFM_surv ~ normal(0, sGeno_surv);
  bGenoNMT_surv ~ normal(0, sGeno_surv);
  bGenoLIJ_surv ~ normal(0, sGeno_surv);
  bGenoRCK_surv ~ normal(0, sGeno_surv);
  
  bSireDamCUR_surv ~ normal(0, sSireDam_surv);
  bSireDamWFM_surv ~ normal(0, sSireDam_surv);
  bSireDamNMT_surv ~ normal(0, sSireDam_surv);
  bSireDamLIJ_surv ~ normal(0, sSireDam_surv);
  bSireDamRCK_surv ~ normal(0, sSireDam_surv);")
    
  }
  
  if (poph2 & nonadd) {
    
    model_priors_surv %<>% 
      stringr::str_c("

  sGenoCUR_surv ~ cauchy(0, 0.1);
  sGenoWFM_surv ~ cauchy(0, 0.1);
  sGenoNMT_surv ~ cauchy(0, 0.1);
  sGenoLIJ_surv ~ cauchy(0, 0.1);
  sGenoRCK_surv ~ cauchy(0, 0.1);
  
  sSireDamCUR_surv ~ cauchy(0, 0.1);
  sSireDamWFM_surv ~ cauchy(0, 0.1);
  sSireDamNMT_surv ~ cauchy(0, 0.1);
  sSireDamLIJ_surv ~ cauchy(0, 0.1);
  sSireDamRCK_surv ~ cauchy(0, 0.1);
  
  bGenoCUR_surv ~ normal(0, sGenoCUR_surv);
  bGenoWFM_surv ~ normal(0, sGenoWFM_surv);
  bGenoNMT_surv ~ normal(0, sGenoNMT_surv);
  bGenoLIJ_surv ~ normal(0, sGenoLIJ_surv);
  bGenoRCK_surv ~ normal(0, sGenoRCK_surv);
  
  bSireDamCUR_surv ~ normal(0, sSireDamCUR_surv);
  bSireDamWFM_surv ~ normal(0, sSireDamWFM_surv);
  bSireDamNMT_surv ~ normal(0, sSireDamNMT_surv);
  bSireDamLIJ_surv ~ normal(0, sSireDamLIJ_surv);
  bSireDamRCK_surv ~ normal(0, sSireDamRCK_surv);")
    
  }
 
  # if (interaction) {
  #   
  #   model_priors_surv %<>% 
  #     stringr::str_replace_all("sGeno", "tauGeno") %>% 
  #     stringr::str_replace_all("sSireDam", "tauSireDam") %>% 
  #     stringr::str_c("\n\tOmegaGeno ~ lkj_corr(2);")
  #   
  #   model_priors_surv %<>% 
  #     stringr::str_replace_all(
  #       "([\\s]*)bGeno([A-Z]{3})_surv([\\s]*~[\\s]*[a-z_]+\\([^\\)]+\\);)", 
  #       "\\1bGeno\\2_surv ~ multi_normal(Zero, quad_form_diag(OmegaGeno, tau\\2_surv));"
  #     )
  #   
  # }  
  # 
  # if (nonadd & interaction) {
  # 
  #   model_priors_surv %<>% 
  #     stringr::str_replace_all(
  #       "([\\s]*)sSireDam([A-Z]{3})_surv([\\s]*~[\\s]*[a-z_]+\\([^\\)]+\\);)",
  #       "\\1tauSireDam\\2_surv\\3"
  #     )
  #   
  #   model_priors_surv %<>% stringr::str_c("\n\tOmegaSireDam_surv ~ lkj_corr(2);")
  #   
  #   model_priors_surv %<>% 
  #     stringr::str_replace_all(
  #       "([\\s]*)bSireDam([A-Z]{3})_surv([\\s]*~[\\s]*[a-z_]+\\([^\\)]+\\);)", 
  #       "\\1bSireDam\\2_surv ~ multi_normal(Zero, quad_form_diag(OmegaSireDam_surv, tauSireDam\\2_surv));"
  #     )
  #   
  # }  
  
  if (interaction) {
    
    model_priors_surv %<>% 
      stringr::str_replace_all("sGeno", "tauGeno") %>% 
      stringr::str_replace_all("sSireDam", "tauSireDam") %>% 
      stringr::str_replace_all(
        "bGeno([A-Z]{3})_surv ~ normal\\(0, tauGeno([A-Z]{0,3})_surv\\);", 
        "bGeno\\1_surv ~ multi_normal(Zero, quad_form_diag(OmegaGeno\\2, tauGeno\\2_surv));"
      ) %>% 
      stringr::str_replace_all(
        "bSireDam([A-Z]{3})_surv ~ normal\\(0, tauSireDam([A-Z]{0,3})_surv\\);", 
        "bSireDam\\1_surv ~ multi_normal(Zero, quad_form_diag(OmegaSireDam\\2, tauSireDam\\2_surv));"
      ) %>% 
      stringr::str_c(
        stringr::str_c(glue::glue("\n\tOmegaGeno{pop} ~ lkj_corr(2);\n", 
                                pop = if (poph2) {pop_levels()} else {""}), 
                       collapse = "\n")) %>%
      stringr::str_c(
        if (nonadd) {
          stringr::str_c(
          glue::glue("\n\tOmegaSireDam{pop} ~ lkj_corr(2);\n", 
                     pop = if (poph2) {pop_levels()} else {""}), collapse = "\n"
          )
          } else {""}
        )
    
  }
  
  model_priors <- stringr::str_c(model_priors_germ, model_priors_surv, "\n")

  ## Model Block (germination target) ----
  
  model_block <- stringr::str_c(model_priors, "\n\t// Target (germination)\n")
  
  for (i in sow_dates$sub_cohort) {
    for (j in levels(data$Pop)) {
      model_block %<>% stringr::str_c(
        make_germ_model(
          pop = j, 
          sub_cohort = i, 
          model = model, 
          nonadd = nonadd
        )
      )
    }
  }
  
  ## Model Block (survival target) ----
  model_block %<>% stringr::str_c("\n\t// Target (survival)\n")
    
  for (pop in levels(data$Pop)) {
    model_block %<>% stringr::str_c(
      make_surv_model(pop = pop, nonadd = nonadd, interaction = interaction)
    )
  }
  
  model_block %<>% stringr::str_c("\n}")
  
  ## Generated Quantities Block (germination) ----
  
  generated_quantities_block_germ <- " generated quantities {\n\n\t// Germination\n\n"
  
  for (i in sow_dates$sub_cohort) {
    for (j in levels(data$Pop)) {
      suffix <- stringr::str_c("_", i, "_", j)
      generated_quantities_block_germ %<>% 
        stringr::str_c("\tvector[n", suffix, "] mu", suffix, "_germ;\n")
    }
  }
  
  generated_quantities_block_germ %<>% stringr::str_c(
    expand.grid("n", sow_dates$sub_cohort, levels(data$Pop)) %>%
      apply(1, stringr::str_c, collapse = "_") %>%
      stringr::str_c(collapse = " + ") %>%
      stringr::str_c(
        "\n\tvector[", ., "] log_lik_germ;\n\tvector[", ., "] predict_germ;"
      )
  )
  
  start <- ""
  
  for (i in sow_dates$sub_cohort) {
    for (j in levels(data$Pop)) {
      generated_quantities_block_germ %<>% 
        stringr::str_c(make_germ_generated_quantities(
          pop = j, 
          sub_cohort = i, 
          start = start,
          model = model, 
          nonadd = nonadd
        ))
      start %<>% stringr::str_c("n", "_", i, "_", j, " + ")
    }
  }
  
  generated_quantities_block_surv <- "\n\n\t// Survival\n\n"
  
  for (pop in levels(data$Pop)) {
    generated_quantities_block_surv %<>% 
      stringr::str_c("\tvector[n", pop, "] theta", pop, ";\n")
  }
  
  generated_quantities_block_surv %<>% stringr::str_c(
    stringr::str_c("n", levels(data$Pop)) %>%
      stringr::str_c(collapse = " + ") %>%
      stringr::str_c("\n\tvector[", ., "] log_lik_surv;\n\tvector[", ., "] predict_surv;"))
  
  start <- ""
  
  for (pop in levels(data$Pop)) {
    generated_quantities_block_surv %<>% 
      stringr::str_c(make_surv_generated_quantities(pop = pop, start = start, nonadd = nonadd, 
                                                    interaction = interaction))
    start %<>% stringr::str_c("n", pop, " + ")
  }
  
  generated_quantities_block <- stringr::str_c(generated_quantities_block_germ,
                                               "\n\n",
                                               generated_quantities_block_surv, 
                                               "\n}\n")
  
  ## Concatenate and Write -----
  
  readr::write_file(
    stringr::str_c(functions_block, data_block, parameters_block, 
                   model_block, generated_quantities_block), 
    file)
  
}


##### Utility functions for writing Stan models ################################

# Add strings for additional parameters depending on model
add_parameters <- function(model, pop) {
  if (model == "discrete_loglogistic") parameters <- 
      ifelse(is.na(pop),
             ", beta, ntd",
             glue::glue(", beta, ntd[{i}]", i = which(pop == pop_levels()))
      )
  if (model == "discrete_lognormal") parameters <- 
      ifelse(is.na(pop),
             ", sigma, ntd",
             glue::glue(", sigma, ntd[{i}]", i = which(pop == pop_levels()))
      )
  parameters
}

add_parameter_types <- function(model) {
  if (model == "discrete_loglogistic") type <- ", real beta, real ntd"
  if (model == "discrete_lognormal") type <- ", real sigma, real ntd"
  type
}

add_parameter_definitions <- function(model) {
  if (model == "discrete_loglogistic") def <- "\t\treal<lower=0> beta;"
  if (model == "discrete_lognormal") def <- "\t\treal<lower=0> sigma;"
  def
}

add_parameter_priors <- function(model) {
  if (model == "discrete_loglogistic") prior <- "\tbeta ~ cauchy(0, 1);\n"
  if (model == "discrete_lognormal") prior <- "\tsigma ~ cauchy(0, 1);\n"
  prior
}

# Suffix for labeling various things in Stan
make_suffix <- function(pop, sub_cohort) {
  suffix <- stringr::str_c("_", sub_cohort, "_", pop)
  suffix
}

# Check model argument
check_model <- function(model) {
  model %<>% match.arg(c("discrete_loglogistic", "discrete_lognormal"))
  model
}

# Check population argument
check_pop <- function(pop) {
  stopifnot(pop %in% pop_levels())
}

pop_levels <- function() {
  c("CUR", "WFM", "NMT", "LIJ", "RCK")
}

##### Miscellaneous ############################################################

# Copy large object and compress it
copy_large_object <- function(x, import = FALSE) {
  
  if (import) {
    
    x %>%
      glue::glue("Importing '{x}'...", x = .) %>%
      crayon::green() %>%
      message()
    
    object_path <- stringr::str_c("r/objects/", x)
    file.copy(from = stringr::str_c("large-objects/", x),
              to = object_path, overwrite = TRUE)
    
    object <- readr::read_rds(object_path)
    readr::write_rds(object, object_path)
    invisible()
    
    
  } else {
    
    x %>%
      glue::glue("Exporting '{x}'...", x = .) %>%
      green() %>%
      message()
    
    object_path <- stringr::str_c("large-objects/", x)
    file.copy(from = stringr::str_c("r/objects/", x),
              to = object_path, overwrite = TRUE)
    
    object <- readr::read_rds(object_path)
    readr::write_rds(object, object_path, compress = "gz")
    invisible()
    
  }
  
}

# Logit and inverse logit functions
logit <- function(x) qlogis(x)
inv_logit <- function(x) plogis(x)

# Function for helping plot hertability
convertY <- function(y, buffer, scal) {
  
  y %<>% 
    multiply_by(1 - 2 * buffer) %>% 
    divide_by(scal) %>% 
    add(buffer) %>%
    grconvertY(from = "npc", "user")
  y
  
}

# Carry forward object for downstream use
carry_forward <- function(object, what, where) {
  
  readr::write_rds(object, dewindowsify_path(stringr::str_c(where, "/", what, ".rds")))
  
}

# Bring forward object from previous analysis
bring_forward <- function(what, where) {
  
  object <- readr::read_rds(dewindowsify_path(stringr::str_c(where, "/", what, ".rds")))
  object 
  
}

# Determine iterations and thinning for Stan
make_stan_options <- function(nsamples, ncores, round = 0) {
  
  thin <- 2 ^ round
  iter <- 2 * nsamples * thin / ncores
  list(thin = thin, iter = iter)
  
} 

# Make labels for pp_check plots
make_pp_subtitle <- function(object_name) {
  
  glue::glue(
    "{distribution}{h2}{Vna}",
    distribution = "Discrete LogNormal",
    h2 = ifelse(stringr::str_detect(object_name, "_commonh2_"),
                "", "\nPopulation-specific heritability"),
    Vna = ifelse(stringr::str_detect(object_name, "_nonadd$"),
                 "\nNonadditive variance", "")
  )
  
}

# De-Windowsify path
dewindowsify_path <- function(path) {
  
  path %<>% stringr::str_replace_all('\\\\', '/')
  path
  
}

# Get population labels from parameter abbreviations
get_labels <- function(string) {
  
  if (stringr::str_detect(string, "CUR")) return("Sweetwater River")
  if (stringr::str_detect(string, "WFM")) return("West Fork Mojave River")
  if (stringr::str_detect(string, "NMT")) return("North Fork Middle Tule River")
  if (stringr::str_detect(string, "LIJ")) return("Little Jamison Creek")
  if (stringr::str_detect(string, "RCK")) return("Rock Creek")
  
}

# Place line break between labels for plotting
place_line_break <- function(s) {
  
  spaces <- stringr::str_locate_all(s, " ")[[1]]
  if (nrow(spaces) == 0L) return(s)
  if (nrow(spaces) == 1L) s %<>% stringr::str_replace(" ", "\n")
  if (nrow(spaces) > 1L) {
    spaces <- spaces[which.min(abs(apply(spaces, 1, mean) - nchar(s) / 2)), ]
    s <- stringr::str_c(substr(s, 1L, spaces["start"] - 1L), "\n", 
                        substr(s, spaces["end"] + 1L, nchar(s)))
  }
  s
  
}

# Get color from palette with custom alpha
get_col <- function(color, alpha) {
  
  x <- col2rgb(color, alpha = TRUE)
  ret <- rgb(x["red", ], x["green", ], x["blue", ], alpha * 255, 
             maxColorValue = 255)
  ret
  
}

# get average emergence time of parents and 'f1'
get_pheno <- function(id, data) {
  
  tmp <- id %>% 
    stringr::str_split("_") %>% 
    extract2(1)
  pop <- tmp[1]
  dam <- stringr::str_c(tmp[1], tmp[2])
  sire <- stringr::str_c(tmp[1], tmp[3])
  rm(tmp)
  pheno_sire <- data %>% 
    dplyr::filter(Sire == sire | Dam == sire) %>% 
    use_series(AvgGerm) %>% 
    mean(na.rm = TRUE)
  pheno_dam <- data %>% 
    dplyr::filter(Sire == dam | Dam == dam) %>% 
    use_series(AvgGerm) %>% 
    mean(na.rm = TRUE)
  pheno_f1 <- data %>% 
    dplyr::filter(Sire == sire & Dam == dam) %>% 
    use_series(AvgGerm) %>% 
    mean(na.rm = TRUE)
  ret <- data.frame(
    id = id, 
    pop = pop, 
    dam = dam, 
    sire = sire, 
    pheno_dam = pheno_dam, 
    pheno_sire = pheno_sire,
    pheno_f1 = pheno_f1
  )
  ret
}

##### Utility functions for manuscript #########################################

# Export objects to ms
export2ms <- function(x, path_export = getOption("path_export", "/export")) {
  
  for (i in 1:length(x)) {
    stopifnot(is.character(x[i]))
    tmp <- eval(parse(text = x[i]))
    path <- normalizePath(stringr::str_c("ms/", path_export)) %>%
      dewindowsify_path()
    if (!dir.exists(path)) dir.create(path)
    readr::write_rds(tmp, stringr::str_c(path, "/", x[i], ".rds"))
  }
  
}

# Import objects to ms
import2ms <- function(path_export = getOption("path_export", "/export")) {
  
  path <- stringr::str_c("ms/", path_export)
  files <- list.files(path, ".rds$")
  object_names <- stringr::str_replace(files, ".rds$", "")
  eval(parse(text = stringr::str_c(
    object_names, " <- readr::read_rds('", path, "/", files, "')")
  ),
  envir = .GlobalEnv)
  
}

# Make statistical significance asterisks
sigStar <- function(pvalue)
{
  if (pvalue >= 0.05) return("n.s.")
  if (pvalue < 0.05 & pvalue >= 0.01) return("*")
  if (pvalue < 0.01 & pvalue >= 0.001) return("**")
  if (pvalue < 0.001) return("***")
}

# Determine a number's order of magnitude
oom <- function(x)
{
  # x should be a vector of numbers
  floor(log10(abs(x)))
}

# Determine number of significant digits to use for rounding in tables
sigDig <- function(x)
{
  # x should be a vector of numeric elements to be rounded to same significant digit
  ret <- round(x, -min(oom(x)))
  return(ret)
}

# Text string for p-value in tables
pval2latex <- function(x)
{
  stem <- x * 10 ^ -oom(x)
  ret <- ifelse(oom(x) > -3, 
                sprintf("%.*f", 3, x), 
                sprintf("%.*f $\\\\times10^{%s}$", 1, stem, oom(x)))
  return(ret)
}



