source("r/header.r")

# Load model
pars <- c("bPop_surv", "bCohortNorth_surv", "tauGeno_surv", "sBlock_surv")
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
  tidyr::pivot_wider(names_from = "parameter",values_from = "value") %>%
  dplyr::ungroup(.draw)

fit_surv_qg <- fit_surv %>%
  dplyr::select(.draw, `tauGeno_surv[1]`, `tauGeno_surv[2]`, 
                dplyr::starts_with("bPop_surv"), 
                dplyr::starts_with("bCohortNorth_surv"), sBlock_surv) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    vpop_south = 4/5 * var(c(`bPop_surv[1]`, `bPop_surv[2]`, `bPop_surv[3]`, 
                             `bPop_surv[4]`, `bPop_surv[5]`)),
    vpop_north = 4/5 * var(c(`bPop_surv[1]` + `bCohortNorth_surv[1]`,
                             `bPop_surv[2]` + `bCohortNorth_surv[2]`,
                             `bPop_surv[3]` + `bCohortNorth_surv[3]`, 
                             `bPop_surv[4]` + `bCohortNorth_surv[4]`, 
                             `bPop_surv[5]` + `bCohortNorth_surv[5]`)),
    mpop_south = mean(c(`bPop_surv[1]`, `bPop_surv[2]`, `bPop_surv[3]`, 
                        `bPop_surv[4]`, `bPop_surv[5]`)),
    mpop_north = mean(c(`bPop_surv[1]` + `bCohortNorth_surv[1]`,
                        `bPop_surv[2]` + `bCohortNorth_surv[2]`,
                        `bPop_surv[3]` + `bCohortNorth_surv[3]`, 
                        `bPop_surv[4]` + `bCohortNorth_surv[4]`, 
                        `bPop_surv[5]` + `bCohortNorth_surv[5]`)),
    va_south = 4 * `tauGeno_surv[1]` ^ 2,
    va_north = 4 * `tauGeno_surv[2]` ^ 2,
    vblock = sBlock_surv ^ 2
  ) %>%
  dplyr::full_join(df_bPop, by = ".draw") %>%
  dplyr::full_join(df_bCohortNorth, by = c(".draw", "pop")) %>%
  dplyr::select(
    bPop_surv, # ..1
    bCohortNorth_surv, # ..2
    va_south, # ..3
    va_north, # ..4
    vblock, # ..5
    mpop_south, # ..6
    mpop_north, # ..7
    vpop_south, # ..8
    vpop_north, # ..9
    pop, # ..10
    .draw # ..11
  ) %>%
  dplyr::filter(.draw < 100) %>%
  purrr::pmap_dfr(~ {
    
    qg_south <- QGglmm::QGparams(mu = ..1, var.a = ..3, var.p = ..3 + ..5, 
                                 model = "binom1.logit", verbose = FALSE)
    qg_block_south <- QGglmm::QGparams(mu = ..1, var.a = ..5,
                                       var.p = ..3 + ..5,
                                       model = "binom1.logit", verbose = FALSE)
    qg_south$var.block.obs <- qg_block_south$var.a.obs
    qg_pop_south <- QGglmm::QGparams(mu = ..6, var.a = ..8,
                                       var.p = ..8 + ..5,
                                       model = "binom1.logit", verbose = FALSE)
    qg_south$var.pop.obs <- qg_pop_south$var.a.obs
    colnames(qg_south) %<>% stringr::str_c(".south")
    
    qg_north <- QGglmm::QGparams(mu = ..1 + ..2, var.a = ..4, var.p = ..4 + ..5,
                                 model = "binom1.logit", verbose = FALSE)
    qg_block_north <- QGglmm::QGparams(mu = ..1 + ..2, var.a = ..4, var.p = ..4 + ..5,
                                       model = "binom1.logit", verbose = FALSE)
    qg_north$var.block.obs <- qg_block_north$var.a.obs
    qg_pop_north <- QGglmm::QGparams(mu = ..7, var.a = ..9,
                                     var.p = ..9 + ..5,
                                     model = "binom1.logit", verbose = FALSE)
    qg_north$var.pop.obs <- qg_pop_north$var.a.obs
    colnames(qg_north) %<>% stringr::str_c(".north")
    
    qg <- dplyr::bind_cols(qg_south, qg_north)
    
    qg$bPop_surv <- ..1
    qg$bCohortNorth_surv <- ..2
    qg$va_south <- ..3
    qg$va_north <- ..4
    qg$vp_south <- ..3 + ..5
    qg$vp_north <- ..3 + ..5
    qg$vpop_south <- ..8
    qg$vpop_north <- ..9
    qg$pop <- ..10
    qg$.draw <- ..11
    qg
    
  })

# Table summarizing variance components
vc_table_surv <- fit_surv_qg %>%
  dplyr::mutate(
    vres_south = var.obs.south - var.a.obs.south - var.block.obs.south,
    vres_north = var.obs.north - var.a.obs.north - var.block.obs.north,
  ) %>%
  dplyr::select(
    .draw, 
    pop,
    va_south = var.a.obs.south,
    va_north = var.a.obs.north,
    vblock_south = var.block.obs.south,
    vblock_north = var.block.obs.north,
    vpop_south = var.pop.obs.south,
    vpop_north = var.pop.obs.north,
    vres_south,
    vres_north,
    h2_south = h2.obs.south,
    h2_north = h2.obs.north
  ) %T>%
  assign(x = "diff_vpop_va", envir = .GlobalEnv) %>%
  tidyr::pivot_longer(-c(pop, .draw), names_to = "parameter") %>%
  dplyr::mutate(
    garden = stringr::str_replace(parameter, 
                                  "^([:alnum:]+)_(south|north)$", "\\2"),
    parameter = stringr::str_replace(parameter, 
                                  "^([:alnum:]+)_(south|north)$", "\\1"),
  ) %>%
  assign("df", ., pos = 1) %>%
  dplyr::group_by(pop, garden, parameter) %>%
  tidybayes::point_interval(.point = median, .interval = tidybayes::qi) %>%
  dplyr::mutate(dplyr::across(value:.upper, signif, digits = 3)) %>%
  dplyr::mutate(
    Population = factor(pop, levels = pop_levels()),
    Garden = stringr::str_to_sentence(garden),
    `Median (95% CI)` = glue::glue("{value} ({.lower}--{.upper})"),
    Parameter = factor(parameter, levels = c("va", "vblock", "vres", "h2"))
  ) %>%
  dplyr::mutate(
    `Median (95% CI)` = stringr::str_replace_all(`Median (95% CI)`, 
                                                  "([0-9].[0-9]{2})e-0([0-9])",
                                                  "$\\1 \\times 10^{-\\2}$")
  ) %>%
  dplyr::arrange(Parameter, Population, Garden) %>%
  dplyr::mutate(
    Parameter = dplyr::case_when(
      parameter == "vpop" ~ "$V_\\text{pop}$",
      parameter == "va" ~ "$V_\\text{A}$",
      parameter == "vblock" ~ "Block",
      parameter == "vres" ~ "$V_\\text{E}$",
      parameter == "h2" ~ "$h^2$"
    )
  ) %>%
  dplyr::select(Population, Garden, Parameter, `Median (95% CI)`) %>%
  dplyr::arrange(Garden, Population, Parameter) %>%
  dplyr::mutate(Population = dplyr::case_when(
    Population == "CUR" ~ "Sweetwater River",
    Population == "WFM" ~ "West Fork Mojave River",
    Population == "NMT" ~ "North Fork Middle Tule River",
    Population == "LIJ" ~ "Little Jamison Creek",
    Population == "RCK" ~ "Rock Creek"
  ))

# Difference between V_pop and V_A for ms
# diff_vpop_va_surv <- diff_vpop_va %>%
#   dplyr::transmute(diff_vpop_va = vpop - va) %>%
#   tidybayes::point_interval(.point = median, .interval = tidybayes::qi) 

# Figure summarizing variance components
df %<>% 
  dplyr::mutate(
    Garden = factor(stringr::str_to_sentence(garden), 
                    levels = c("South", "North")),
    par = factor(dplyr::case_when(
      parameter == "vpop" ~ "italic(V)[pop]",
      parameter == "va" ~ "italic(V)[A]",
      parameter == "vblock" ~ "Block", 
      parameter == "vres" ~ "italic(V)[E]", 
      parameter == "h2" ~ "italic(h) ^ 2"
    ), levels = c("Block", "italic(V)[E]", "italic(V)[A]",
                  "italic(h) ^ 2", "italic(V)[pop]"))
  )

df_summary <- df %>%
  dplyr::group_by(pop, Garden, par) %>%
  tidybayes::point_interval(value, .width = c(0.8, 0.95), .point = median,
                            .interval = tidybayes::qi) %>%
  dplyr::filter(!(pop != "CUR" & par == "italic(V)[pop]")) %>%
  dplyr::mutate(pop = factor(pop, levels = c(pop_levels(), "all")))

df_summary[df_summary$par == "italic(V)[pop]", "pop"] <- "all"

ggplot(df, aes(pop, value, color = pop)) +
  facet_grid(Garden ~ par, scales = "free_x", space = "free_x", 
             as.table = FALSE, labeller = label_parsed) +
  # geom_violin(trim = FALSE, scale = "width", adjust = 2) +
  geom_linerange(
    data = dplyr::filter(df_summary, .width == 0.95), 
    mapping = aes(ymin = .lower, ymax = .upper), 
    position = position_dodge(width = 0.9)
  ) +
  geom_linerange(
    data = dplyr::filter(df_summary, .width == 0.8), 
    mapping = aes(ymin = .lower, ymax = .upper), 
    position = position_dodge(width = 0.9), size = 1.2
  ) +
  geom_point(
    data = dplyr::filter(df_summary, .width == 0.8), 
    position = position_dodge(width = 0.9), size = 3
  ) +
  scale_color_manual(
    values = c(palette(), "grey50"),
    labels = pop_levels() %>% sapply(get_labels) %>% sapply(place_line_break) %>%
      c(all = "Among\npopulation"),
    name = "Population\nof origin"
  ) +
  scale_y_log10(labels = scales::label_math(expr = 10^.x, format = log10) ) +
  xlab("Population") +
  ylab(expression(paste("Variance or Heritability (", log[10], "-scale)"))) +
  ggtitle("Winter survival") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 14),
    legend.position = "bottom",
    legend.text = element_text(size = 12, hjust = 0.5),
    legend.title = element_text(size = 14, hjust = 0.5),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12)
  )

ggsave("ms/figures/h2-surv.pdf", width = 7.5, height = 7.5, units = "in")

export2ms("vc_table_surv")
