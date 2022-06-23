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

bPop = dplyr::full_join(df_bPop, df_bCohortNorth, by = c(".draw", "pop")) |>
  dplyr::mutate(north = bPop_surv + bCohortNorth_surv, south = bPop_surv) |>
  dplyr::select(.draw, pop, south, north) |>
  tidyr::pivot_longer(cols = c("south", "north"), names_to = "garden") |>
  dplyr::rename(pop1 = pop, garden1 = garden, value1 = value)

letter_table <- bPop |>
  dplyr::select(tidyr::starts_with(".")) |>
  dplyr::distinct() |>
  tidyr::crossing(pop2 = pop_levels(), garden2 = c("south", "north")) |>
  dplyr::left_join(bPop) |>
  dplyr::left_join(dplyr::rename(bPop, pop2 = pop1, garden2 = garden1, value2 = value1)) |>
  dplyr::mutate(diff_surv = value1 - value2) |>
  dplyr::group_by(pop1, garden1, pop2, garden2) |>
  ggdist::point_interval(value1, diff_surv, .point = median, .interval = tidybayes::qi) |>
  dplyr::mutate(sig_diff = sign(diff_surv.lower) == sign(diff_surv.upper)) |>
  dplyr::filter(diff_surv != 0) |>
  dplyr::select(pop1, garden1, pop2, garden2, value1, sig_diff) |>
  dplyr::arrange(value1)
  
df_letters <- tidyr::crossing(
  pop = pop_levels(),
  garden = c("south", "north")
) |>
  dplyr::mutate(
    pop = factor(pop, levels = pop_levels()),
    garden = factor(garden, levels = c("south", "north"))
  ) |>
  dplyr::arrange(pop, garden) |>
  dplyr::mutate(
    letter = c("f", "cd", "ef", "c", "ef", "d", "e", "b", "cd", "a"),
    p_surv = c(1.05, 1, 1.05, 0.975, 1.05, 1, 1.025, 0.975, 1, 0.85)
  )

# print(letter_table, n = 81)
# 1 CUR   south f
# 2 CUR   north cd
# 3 WFM   south ef
# 4 WFM   north c
# 5 NMT   south ef
# 6 NMT   north d
# 7 LIJ   south e
# 8 LIJ   north b
# 9 RCK   south cd
# 10 RCK   north a

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
  geom_text(data = dplyr::filter(df_letters, garden == "south"), aes(label = letter), nudge_x = -0.25) +
  geom_text(data = dplyr::filter(df_letters, garden == "north"), aes(label = letter), nudge_x = 0.25) +
  scale_x_discrete(
    labels = pop_levels() %>% sapply(get_labels) %>% sapply(place_line_break)
  ) +
  scale_fill_manual(values = palette(), guide = "none") +
  scale_linetype(name = "Garden") +
  xlab("Source Population") +
  ylab("Winter Survival") +
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

gp <- plot_grid(mean_germ, mean_surv, ncol = 1, labels = "auto")
ggsave("ms/figures/mean-traits.pdf", plot = gp, width = 5, height = 10,
       units = "in")
