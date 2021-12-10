source("r/header.R")

# Load model and data -----
fit <- readr::read_rds("r/objects/fit.rds")
data <- readr::read_rds("processed-data/germination.rds")

# Extract breeding values -----

germ_pop <- fit$draws("bPop_germ") %>%
  posterior::as_draws_df() %>%
  dplyr::select(dplyr::starts_with("bPop_germ"), .draw) %>%
  tidyr::pivot_longer(-.draw) %>%
  dplyr::mutate(
    pop = stringr::str_extract(name, "[1-5]{1}"),
    pop = pop_levels()[as.numeric(pop)]
  ) %>%
  dplyr::rename(bPop_germ = value) %>%
  dplyr::select(-name)

germ_ind <- fit$draws(glue::glue("bGeno{pop}_germ", pop = pop_levels())) %>%
  posterior::as_draws_df() %>%
  dplyr::select(dplyr::starts_with("bGeno"), .draw) %>%
  tidyr::pivot_longer(-.draw) %>%
  dplyr::mutate(
    pop = stringr::str_extract(name, "[A-Z]{3}"),
    ind = stringr::str_extract(name, "[0-9]{1,2}")
  ) %>%
  dplyr::rename(bGeno_germ = value) %>%
  dplyr::select(-name) %>%
  dplyr::full_join(germ_pop, by = c(".draw", "pop")) %>%
  dplyr::mutate(DaysToGerm = 4 + exp(bPop_germ + bGeno_germ))

surv_pop_south <- fit$draws("bPop_surv") %>%
  posterior::as_draws_df() %>%
  dplyr::select(dplyr::starts_with("bPop_surv"), .draw) %>%
  tidyr::pivot_longer(-.draw) %>%
  dplyr::mutate(
    pop = stringr::str_extract(name, "[1-5]{1}"),
    pop = pop_levels()[as.numeric(pop)]
  ) %>%
  dplyr::select(.draw, pop, bPop_surv = value)

df_bCohortNorth <- fit$draws("bCohortNorth_surv") %>%
  posterior::as_draws_df() %>%
  dplyr::select(.draw, dplyr::starts_with("bCohortNorth")) %>%
  tidyr::pivot_longer(-.draw, names_to = "parameter") %>%
  dplyr::mutate(
    pop = stringr::str_extract(parameter, "[0-9]{1,2}"),
    pop = pop_levels()[as.numeric(pop)],
    pop = factor(pop, pop_levels()),
    parameter = "bCohortNorth_surv"
  ) %>%
  dplyr::group_by(.draw) %>%
  tidyr::pivot_wider(names_from = "parameter", values_from = "value") %>%
  dplyr::ungroup(.draw)

surv_pop_north <- surv_pop_south %>%
  dplyr::full_join(df_bCohortNorth, by = c(".draw", "pop")) %>%
  dplyr::transmute(.draw = .draw, pop = pop, 
                   bPop_surv = bPop_surv + bCohortNorth_surv)
  
surv_pop <- dplyr::full_join(
  dplyr::rename(surv_pop_south, south = bPop_surv),
  dplyr::rename(surv_pop_north, north = bPop_surv),
  by = c(".draw", "pop")
) %>%
  tidyr::pivot_longer(-c(.draw, pop), names_to = "garden", 
                      values_to = "logit_p_surv")

surv_ind <- fit$draws(glue::glue("bGeno{pop}_surv", pop = pop_levels())) %>%
  posterior::as_draws_df() %>%
  dplyr::select(.draw, dplyr::starts_with("bGeno")) %>%
  tidyr::pivot_longer(-.draw) %>%
  dplyr::mutate(
    pop = stringr::str_extract(name, "[A-Z]{3}"),
    ind = stringr::str_replace(name, "bGeno[A-Z]{3}_surv\\[([0-9]{1,2}),([1-2]{1})\\]$", "\\1"),
    garden = stringr::str_replace(name, "bGeno[A-Z]{3}_surv\\[([0-9]{1,2}),([1-2]{1})\\]$", "\\2"),
    garden = dplyr::case_when(
      garden == "1" ~ "south",
      garden == "2" ~ "north"
    )
  ) %>%
  dplyr::full_join(surv_pop, by = c(".draw", "pop", "garden")) %>%
  dplyr::select(.draw, garden, pop, ind, bGeno_surv = value, 
                bPop_surv = logit_p_surv, -name)

# Join Germination and Survival posteriors ----
pop <- germ_pop %>%
  dplyr::mutate(DaysToGerm = 4 + exp(bPop_germ)) %>%
  dplyr::full_join(surv_pop, by  = c(".draw", "pop"))

ind <- dplyr::full_join(germ_ind, surv_ind, by = c(".draw", "pop", "ind"))

# Linear Regressions ----

df_pop_lm <- pop %>% 
  dplyr::select(DaysToGerm, logit_p_surv, .draw, pop, garden) %>%
  dplyr::group_by(.draw, garden) %>%
  dplyr::summarise(
    slope = cor(DaysToGerm, logit_p_surv) * sd(logit_p_surv) / sd(DaysToGerm),
    intercept = mean(logit_p_surv) - slope * mean(DaysToGerm)
  )

df_ind_lm <- ind %>% 
  dplyr::select(bGeno_germ, bGeno_surv, .draw, pop, ind, garden) %>%
  dplyr::group_by(.draw, garden) %>%
  dplyr::summarise(
    slope = cor(bGeno_germ, bGeno_surv) * sd(bGeno_surv) / sd(bGeno_germ),
    intercept = mean(bGeno_surv) - slope * mean(bGeno_germ)
  )

# Quadratic Regressions ----

df_pop_lm2 <- pop %>% 
  dplyr::select(DaysToGerm, logit_p_surv, .draw, pop, garden) %>%
  dplyr::group_by(.draw, garden) %>%
  purrrlyr::by_slice(~{
    fit <- lm(logit_p_surv ~ poly(DaysToGerm, 2, raw = TRUE), data = .x)
    b <- coef(fit)
    tibble::tibble(b0 = b[1], b1 = b[2], b2 = b[3])
  }, .collate = "rows")

df_pop_lm2 %>%
  dplyr::group_by(garden) %>%
  tidybayes::point_interval(b1)
  
df_ind_lm2 <- ind %>% 
  dplyr::select(bGeno_germ, bGeno_surv, .draw, pop, ind, garden) %>%
  dplyr::group_by(.draw, garden) %>%
  purrrlyr::by_slice(~{
    fit <- lm(bGeno_surv ~ poly(bGeno_germ, 2, raw = TRUE), data = .x)
    b <- coef(fit)
    tibble::tibble(b0 = b[1], b1 = b[2], b2 = b[3])
  }, .collate = "rows")

df_ind_lm2 %>%
  dplyr::group_by(garden) %>%
  tidybayes::point_interval(b2)

# Figure ----

df_pop_lines <- tidyr::crossing(
    .draw = 1:max(df_pop_lm$.draw),
    DaysToGerm = seq(min(pop$DaysToGerm), max(pop$DaysToGerm), length.out = 1e2)
  ) %>%
  dplyr::full_join(df_pop_lm, by = ".draw") %>%
  dplyr::mutate(logit_p_surv = intercept + slope * DaysToGerm) %>%
  dplyr::group_by(garden, DaysToGerm) %>%
  tidybayes::point_interval(logit_p_surv, .point = median, 
                            .interval = tidybayes::qi, .width = 0.95) %>%
  dplyr::mutate(garden = factor(garden, levels = c("south", "north")))

df_pop_points <- pop %>%
  dplyr::group_by(pop, garden) %>%
  tidybayes::point_interval(DaysToGerm, logit_p_surv, .point = median, 
                            .interval = tidybayes::qi, .width = 0.95) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    pop = factor(pop, levels = pop_levels()),
    garden = factor(garden, levels = c("south", "north")),
    .lower = logit_p_surv.lower,
    .upper = logit_p_surv.upper
  )

df_ind_points <- ind %>%
  dplyr::mutate(logit_p_surv = bPop_surv + bGeno_surv) %>%
  dplyr::group_by(pop, ind, garden) %>%
  dplyr::summarize(dplyr::across(c(DaysToGerm, logit_p_surv), .fns = median)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    .lower = NA, 
    .upper = NA, 
    pop = factor(pop, levels = pop_levels())
  )
  
ggplot(df_pop_lines, 
       aes(DaysToGerm, inv_logit(logit_p_surv), fill = garden, group = garden,
           ymin = inv_logit(.lower), ymax = inv_logit(.upper))
       ) +
  geom_ribbon() +
  geom_line(size = 1.2, lineend = "round") +
  geom_line(mapping = aes(y = inv_logit(.lower)), color = "grey25") +
  geom_line(mapping = aes(y = inv_logit(.upper)), color = "grey25") +
  geom_point(data = df_ind_points, mapping = aes(color = pop), size = 2,
             alpha = 0.5) +
  geom_linerange(data = df_pop_points) +
  geom_errorbarh(
    data = df_pop_points,
    mapping = aes(xmin = DaysToGerm.lower, xmax = DaysToGerm.upper)
  ) +
  geom_point(data = df_pop_points, mapping = aes(color = pop), size = 3) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(
    values = palette()[1:5], name = "population",
    labels = sapply(pop_levels(), 
                    function(.x) place_line_break(get_labels(.x)))
  ) +
  scale_fill_manual(values = c("grey", "white")) +
  xlab("Days to Germination") +
  ylab("Winter Survival") +
  guides(color = guide_legend(keyheight = 3)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text.align = 0.5
  )

ggsave("ms/figures/selection.pdf", width = 6.5, height = 5, units = "in")


# Table - should be combined with genotypic results ----
df_pop_lm %<>%
  dplyr::group_by(garden) %>%
  tidyr::pivot_longer(c(-.draw, -garden), names_to = "parameter") %>%
  dplyr::group_by(garden, parameter) %>%
  tidybayes::point_interval(
    value, 
    .width = 0.95, 
    .point = median,
    .interval = tidybayes::qi
  ) %>%
  dplyr::mutate(
    dplyr::across(value:.upper, signif, 3),
    Garden = dplyr::case_when(
      garden == "north" ~ "North",
      garden == "south" ~ "South"
    ),
  `Median (95% CI)` = glue::glue("{value} ({.lower} -- {.upper})"),
  ) %>%
  dplyr::select(Garden, Parameter = parameter, `Median (95% CI)`)

df_ind_lm %<>%
  dplyr::group_by(garden) %>%
  tidyr::pivot_longer(c(-.draw, -garden), names_to = "parameter") %>%
  dplyr::group_by(garden, parameter) %>%
  tidybayes::point_interval(
    value, 
    .width = 0.95, 
    .point = median,
    .interval = tidybayes::qi
  ) %>%
  dplyr::mutate(
    dplyr::across(value:.upper, signif, 3),
    Garden = dplyr::case_when(
      garden == "north" ~ "North",
      garden == "south" ~ "South"
    ),
    `Median (95% CI)` = glue::glue("{value} ({.lower} -- {.upper})"),
  ) %>%
  dplyr::select(Garden, Parameter = parameter, `Median (95% CI)`)

export2ms(c("df_pop_lm", "df_ind_lm"))
