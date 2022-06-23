source("r/header.r")

# Extract likelihoods ----
models = list.files(
  path = "r/objects", 
  pattern = "^lognormal_[01]{1}_1.rds"
) %>%
  stringr::str_c("r/objects/", .)

# Germination model ----

likelihoods = models %>%
  parallel::mclapply(function(.x) {
    readr::read_rds(.x)
  }, mc.cores = min(c(parallel::detectCores(), length(.)))) %>%
  parallel::mclapply(function(.x) {.x$draws("log_lik_germ")},
                     mc.cores = min(c(parallel::detectCores(), length(.))))

r_eff = likelihoods %>%
  purrr::map(~ {
    .x = exp(posterior::as_draws_array(.x))
    loo::relative_eff(.x, cores = parallel::detectCores())
  })

names(likelihoods) = names(r_eff) = models

## Compare models using LOOIC -----
loo_table_germ = purrr::map2(likelihoods, r_eff, ~ {
  loo::loo(posterior::as_draws_array(.x), r_eff = .y, 
           cores = parallel::detectCores())
}) %>%
  loo::loo_compare() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("model") %>%
  dplyr::mutate(
    preferred_model = model == "r/objects/lognormal_0_1.rds",
    signif_diff = (abs(elpd_diff) > 2 * se_diff)
  )

best_model = ifelse(
  loo_table_germ$signif_diff[loo_table_germ$preferred_model],
  dplyr::first(loo_table_germ$model),
  loo_table_germ$model[loo_table_germ$preferred_model]
)

# Survival model ----

likelihoods = models %>%
  parallel::mclapply(function(.x) {
    readr::read_rds(.x)
  }, mc.cores = min(c(parallel::detectCores(), length(.)))) %>%
  parallel::mclapply(function(.x) {.x$draws("log_lik_surv")},
                     mc.cores = min(c(parallel::detectCores(), length(.))))

r_eff = likelihoods %>%
  purrr::map(~ {
    .x = exp(posterior::as_draws_array(.x))
    loo::relative_eff(.x, cores = parallel::detectCores())
  })

names(likelihoods) = names(r_eff) = models

## Compare models using LOOIC -----
loo_table_surv = purrr::map2(likelihoods, r_eff, ~ {
  loo::loo(posterior::as_draws_array(.x), r_eff = .y, 
           cores = parallel::detectCores())
}) %>%
  loo::loo_compare() %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("model") %>%
  dplyr::mutate(
    preferred_model = model == "r/objects/lognormal_0_1.rds",
    signif_diff = (abs(elpd_diff) > 2 * se_diff)
  )

best_model = ifelse(
  loo_table_surv$signif_diff[loo_table_surv$preferred_model],
  dplyr::first(loo_table_surv$model),
  loo_table_surv$model[loo_table_surv$preferred_model]
)

# Both analyses agree there is no evidence for population-specific genetic
# variance components. 

file.copy("r/objects/lognormal_0_1.rds", "r/objects/fit.rds", overwrite = TRUE)
cat("r/objects/fit.rds filter=lfs diff=lfs merge=lfs -text\n",
    file = ".gitattributes", append = TRUE)
