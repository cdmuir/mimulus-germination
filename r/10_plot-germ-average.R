source("r/header.r")

# Load model
pars_pop <- glue::glue("bPop_germ[{n}]", n = 1:5)
pars_ind <-  stringr::str_c("bGeno", pop_levels(), "_germ")
fit_germ <- readr::read_rds("r/objects/fit.rds") 

df <- fit_germ$draws(c(pars_pop, pars_ind)) %>% 
  posterior::as_draws_df()

df_pop <- df %>%
  dplyr::select(.draw, pars_pop) %>%
  tidyr::pivot_longer(-.draw, names_to = "parameter") %>%
  dplyr::mutate(
    pop = stringr::str_replace(parameter, "^bPop_germ\\[([1-5])\\]$", "\\1"),
    pop = pop_levels()[as.numeric(pop)],
    pop = factor(pop, pop_levels()),
    parameter = "bPop_germ"
  ) %>%
  dplyr::group_by(.draw) %>%
  tidyr::spread(parameter, value) %>%
  dplyr::ungroup(.draw) 

df_ind <- df %>%
  dplyr::select(.draw, starts_with("bGeno")) %>%
  tidyr::pivot_longer(-.draw, names_to = "parameter") %>%
  dplyr::mutate(
    pop = stringr::str_extract(parameter, "[A-Z]{3}"),
    pop = factor(pop, pop_levels()),
    ind = stringr::str_extract(parameter, "[0-9]{1,2}"),
    parameter = "bGeno"
  ) %>%
  dplyr::group_by(.draw) %>%
  tidyr::pivot_wider(names_from = "parameter", values_from = "value") %>%
  dplyr::ungroup(.draw) %>%
  dplyr::full_join(df_pop, by = c(".draw", "pop")) %>%
  dplyr::mutate(DaysToGerm = 4 + exp(bPop_germ + bGeno)) %>%
  dplyr::group_by(pop, ind) %>%
  tidybayes::point_interval(DaysToGerm, .point = median, .interval = tidybayes::qi)
  
df_pop %<>% dplyr::mutate(DaysToGerm = 4 + exp(bPop_germ))
  
mean_germ <- ggplot(df_pop, aes(pop, DaysToGerm, fill = pop)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 2, draw_quantiles = 0.5,
              show.legend = FALSE) +
  # geom_violin(trim = FALSE, scale = "width", fill = "grey", adjust = 2,
              # draw_quantiles = 0.5) +
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
  #              geom = "tile", size = 1, height = 0.05, width = 0.5,
  #              fill = "grey50") +
  geom_point(data = df_ind, position = position_jitter(width = 0.1, height = 0),
             show.legend = FALSE) + 
  scale_x_discrete(
    labels = pop_levels() %>% sapply(get_labels) %>% sapply(place_line_break)
  ) +
  scale_fill_manual(values = palette()) +
  xlab("Population") +
  ylab("Days to Germination") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 30, vjust = 0.75),
    axis.text = element_text(size = 12),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12)
  )

readr::write_rds(mean_germ, "r/objects/mean_germ.rds")

# ggsave("ms/figures/mean-germ.pdf", width = 5, height = 5, units = "in")
