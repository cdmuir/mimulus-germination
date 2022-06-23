source("r/header.r")

# Posterior predictive checks -----

sow_dates = readr::read_rds("r/objects/sow_dates.rds")

# Data
data <- readr::read_rds("r/objects/germ_stan.rds")

y <- tibble::tibble(
  DaysToGerm = data$DaysToGerm,
  pop = factor(data$Pop_germ, levels = rev(pop_levels())),
  cohort = data$Cohort_germ
) %>%
  dplyr::group_by(pop, cohort, DaysToGerm) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")

# Model predictions
fit <- readr::read_rds("r/objects/fit.rds")

# Not plotted but useful for checking
df_expected <- tidyr::crossing(
  tidyr::nesting(
    pop = factor(pop_levels(), levels = rev(pop_levels())), 
    mu = rev(fit$summary("bPop_germ")$median),
    sigma = rev(fit$summary("sigma")$median)
  ),
  cohort = c(0, 1),
  b_cohort = fit$summary("bCohortNorth_germ")$median,
  DaysToGerm = 5:max(y$DaysToGerm)
) %>%
  dplyr::full_join(
    y %>%
      dplyr::group_by(pop, cohort) %>%
      dplyr::summarize(n_total = sum(n))
  ) %>%
  dplyr::mutate(
    prob_mass = plnorm(DaysToGerm - 4, mu + b_cohort * cohort, sigma) - 
      plnorm(DaysToGerm - 5, mu + b_cohort * cohort, sigma),
    n = n_total * prob_mass
  )

# Posterior predictions
yrep <- fit$draws("predict_germ") %>%
  posterior::as_draws_df() %>%
  dplyr::select(-.chain, -.iteration) %>%
  tidyr::pivot_longer(-.draw, values_to = "DaysToGerm") %>%
  dplyr::mutate(i = as.numeric(stringr::str_extract(name, "[0-9]+"))) %>%
  dplyr::full_join(
    tidyr::crossing(
      sub_cohort = sow_dates$sub_cohort,
      pop = factor(pop_levels(), levels = pop_levels())
    ) %>%
      purrr::pmap_dfr(~{
        n <- data[[glue::glue("n_{sub_cohort}_{pop}", sub_cohort = ..1, pop = ..2)]]
        tibble::tibble(sub_cohort = rep(..1, n), pop = rep(..2, n))
      }) %>%
      dplyr::mutate(
        i = dplyr::row_number(),
        cohort = as.numeric(stringr::str_detect(sub_cohort, "north"))
      ), by = "i"
  ) %>%
  dplyr::group_by(.draw, pop, cohort, DaysToGerm) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(pop, cohort, DaysToGerm) %>%
  tidybayes::point_interval(n, .width = 0.9, .interval = tidybayes::qi)

# Plot

pop_labels = pop_levels() %>% 
  rev() %>% 
  sapply(get_labels) %>% 
  sapply(place_line_break)
cohort_labels = c(`0` = "South cohort", `1` = "North cohort")

gp <- ggplot(y, aes(DaysToGerm, n, fill = "Observations")) +
  facet_grid(rows = vars(pop), cols = vars(cohort),
             labeller = labeller(pop = pop_labels, cohort = cohort_labels)) +
  geom_col(color = "white") +
  geom_pointrange(
    data = yrep, 
    mapping = aes(DaysToGerm, n, ymin = .lower, ymax = .upper, 
                  color = "Posterior\npredictions"), show.legend = TRUE,
    size = 0.75
    ) +
  xlab("Days to Germination") +
  ylab("Count") +
  ggtitle(label = "Posterior Predictive Check") + #, subtitle = subtitle) +
  scale_fill_manual(name = NULL, values = "grey50") +
  scale_color_manual(name = NULL, values = "black") +
  guides(fill = guide_legend(override.aes = list(color = "grey50"))) +
  theme(
    axis.title = element_text(size = 12),
    legend.background = element_rect(fill = NA),
    legend.justification = c(1, 1),
    legend.position = "top",
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )

ggplot2::ggsave("ms/figures/pp_check_germ.pdf", plot = gp, width = 6.5, 
                height = 9, useDingbats = FALSE)
