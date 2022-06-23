source("r/header.R")

fit <- readr::read_rds("r/objects/fit.rds")
fit_summary <- fit$summary()

fit_summary %>%
  dplyr::filter(rhat > 1.05)

fit_summary %>%
  dplyr::filter(ess_tail < 1e3)
