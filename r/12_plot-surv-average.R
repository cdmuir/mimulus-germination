source("r/header.r")

pars <- c("bPop_surv", "bCohortNorth_surv", "tauGeno_surv")
fit_surv <- readr::read_rds("r/objects/fit.rds")$draws(pars) %>%
  posterior::as_draws_df()

# Calculate additive genetic variance posterior
df_bPop <- fit_surv %>%
  dplyr::select(.draw, dplyr::starts_with("bPop")) %>%
  tidyr::pivot_longer(-.draw, names_to = "parameter") %>%
  dplyr::mutate(
    pop = stringr::str_extract(parameter, "[0-9]{1,2}"),
    pop = pop_levels()[as.numeric(pop)],
    pop = factor(pop, pop_levels()),
    parameter = "bPop_surv"
  ) %>%
  dplyr::group_by(.draw) %>%
  tidyr::pivot_wider(names_from = "parameter",values_from = "value") %>%
  dplyr::ungroup(.draw)

df_bCohortNorth <- fit_surv %>%
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

dplyr::full_join(df_bPop, df_bCohortNorth, by = c(".draw", "pop"))

pars <- glue::glue("bGeno{pop}_surv", pop = pop_levels())
fit_surv <- readr::read_rds("r/objects/fit.rds")$draws(pars) %>%
  posterior::as_draws_df()

df_ind <- fit_surv %>%
  dplyr::select(.draw, dplyr::starts_with("bGeno")) %>%
  tidyr::pivot_longer(-.draw) %>%
  dplyr::mutate(
    pop = stringr::str_extract(name, "[A-Z]{3}"),
    geno = stringr::str_replace(name, "bGeno[A-Z]{3}_surv\\[([0-9]{1,2}),([1-2]{1})\\]$", "\\1"),
    garden = stringr::str_replace(name, "bGeno[A-Z]{3}_surv\\[([0-9]{1,2}),([1-2]{1})\\]$", "\\2"),
    garden = dplyr::case_when(
      garden == "1" ~ "south",
      garden == "2" ~ "north"
    )
  ) %>%
  dplyr::full_join(df_bPop, by = c(".draw", "pop")) %>%
  dplyr::full_join(df_bCohortNorth, by = c(".draw", "pop")) %>%
  dplyr::group_by(garden, pop, geno) %>%
  dplyr::mutate(p_surv = plogis(bPop_surv + value + 
                  (garden == "north") * bCohortNorth_surv)) %>% 
  tidybayes::point_interval(p_surv) %>%
  dplyr::mutate(garden = factor(garden, levels = c("south", "north")))

df_pop <- df_bPop %>%
  dplyr::full_join(df_bCohortNorth, by = c(".draw", "pop")) %>%
  tidyr::crossing(garden = c("south", "north")) %>%
  dplyr::mutate(garden = factor(garden, levels = c("south", "north"))) %>%
  dplyr::mutate(p_surv = plogis(bPop_surv + 
                                  (garden == "north") * bCohortNorth_surv))

mean_surv <- ggplot(df_pop, aes(pop, p_surv, linetype = garden, fill = pop)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 2, draw_quantiles = 0.5) +
  geom_point(data = df_ind, 
             position = position_jitterdodge(
               jitter.width = 0.5, 
               jitter.height = 0, 
               dodge.width = 0.9),
             show.legend = FALSE
             ) + 
  scale_x_discrete(
    labels = pop_levels() %>% sapply(get_labels) %>% sapply(place_line_break)
  ) +
  scale_fill_manual(values = palette(), guide = "none") +
  xlab("Population") +
  ylab("Winter Survival") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 30, vjust = 0.75),
    axis.text = element_text(size = 12),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 14),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12)
  )

mean_germ <- readr::read_rds("r/objects/mean_germ.rds")

plot_grid(mean_germ, mean_surv, ncol = 1, labels = "auto")
ggsave("ms/figures/mean-traits.pdf", width = 5, height = 10, units = "in")
