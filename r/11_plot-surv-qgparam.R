source("r/header.r")

# Load model
pars <- c("bPop_surv", "bCohortNorth_surv", "tauGeno_surv", "sBlock_surv")
fit_surv <- readr::read_rds("r/objects/fit.rds")$draws(pars) %>%
  posterior::as_draws_df()

# Calculate additive genetic variance posterior
df_pop <- fit_surv %>%
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

# need to adjust for pop-specific bCohortNorth_surv
fit_surv %>%
  dplyr::select(-dplyr::starts_with("bPop")) %>%
  dplyr::mutate(
    va_south = 4 * `tauGeno_surv[1]` ^ 2,
    va_north = 4 * `tauGeno_surv[2]` ^ 2,
    vblock = sBlock_surv ^ 2
  ) %>%
  dplyr::full_join(df_pop, by = ".draw") %>%
  dplyr::select(bPop_surv, bCohortNorth_surv, va_south, va_north, vblock, pop,
                .draw) %>%
  purrr::pmap_dfr(~ {
    
    qg_south <- QGglmm::QGparams(mu = ..1, var.a = ..3, var.p = ..3 + ..5, 
                                 model = "binom1.logit", verbose = FALSE)
    qg_block_south <- QGglmm::QGparams(mu = ..1, var.a = ..5, 
                                       var.p = ..3 + ..5, 
                                 model = "binom1.logit", verbose = FALSE)
    qg_south$var.block.obs <- qg_block_south$var.a.obs
    colnames(qg_south) %<>% stringr::str_c(".south")
    
    qg_north <- QGglmm::QGparams(mu = ..1 + ..2, var.a = ..4, var.p = ..4 + ..5, 
                                 model = "binom1.logit", verbose = FALSE)
    qg_block_north <- QGglmm::QGparams(mu = ..1 + ..2, var.a = ..4, var.p = ..4 + ..5, 
                                 model = "binom1.logit", verbose = FALSE)
    qg_north$var.block.obs <- qg_block_north$var.a.obs
    colnames(qg_north) %<>% stringr::str_c(".north")
    
    qg <- dplyr::bind_cols(qg_south, qg_north)
    
    qg$bPop_surv <- ..1
    qg$bCohortNorth_surv <- ..2
    qg$va_south <- ..3
    qg$va_north <- ..4
    qg$vp_south <- ..3 + ..5
    qg$vp_north <- ..3 + ..5
    qg$pop <- ..6
    qg$.draw <- ..7
    qg
    
  })

# Table summarizing variance components
vc_table_surv <- fit_surv %>%
  dplyr::mutate(
    vres_south = var.obs.south - var.a.obs.south - var.block.obs.south,
    vres_north = var.obs.north - var.a.obs.north - var.block.obs.north,
  ) %>%
  dplyr::select(
    .iter, 
    pop,
    va_south = var.a.obs.south,
    va_north = var.a.obs.north,
    vblock_south = var.block.obs.south,
    vblock_north = var.block.obs.north,
    vres_south,
    vres_north,
    h2_south = h2.obs.south,
    h2_north = h2.obs.north
  ) %>%
  tidyr::gather(parameter, value, -pop, -.iter) %>%
  dplyr::mutate(
    garden = stringr::str_replace(parameter, 
                                  "^([:alnum:]+)_(south|north)$", "\\2"),
    parameter = stringr::str_replace(parameter, 
                                  "^([:alnum:]+)_(south|north)$", "\\1"),
  ) %>%
  assign("df", ., pos = 1) %>%
  dplyr::group_by(.iter, pop, garden) %>%
  tidyr::spread(parameter, value) %>%
  dplyr::group_by(pop, garden) %>%
  tidybayes::point_interval(va, vblock, vres, h2, 
                            .interval = tidybayes::hdci) %>%
  dplyr::select_at(dplyr::vars(-dplyr::starts_with("."))) %>%
  tidyr::gather(key, value, -pop, -garden) %>%
  dplyr::mutate(
    parameter = stringr::str_replace(key, "^([:alnum:]+).*[a-z]*$", "\\1"),
    type = stringr::str_replace(key, "^[:alnum:]+[.]*([a-z]*)$", "\\1"),
    type = ifelse(nchar(type) == 0L, "point", type)
  ) %>%
  dplyr::select(-key) %>%
  dplyr::group_by(pop, garden) %>%
  tidyr::spread(type, value) %>%
  dplyr::ungroup() %>%
  dplyr::mutate_if(is.numeric, signif, digits = 3) %>%
  dplyr::mutate(
    Population = factor(pop, levels = pop_levels()),
    Garden = stringr::str_to_sentence(garden),
    `Median (95% HDI)` = glue::glue("{point} ({lower}--{upper})",
                                    point = point, lower = lower, upper = upper),
    Parameter = factor(parameter, levels = c("va", "vblock", "vres", "h2"))
  ) %>%
  dplyr::mutate(
    `Median (95% HDI)` = stringr::str_replace_all(`Median (95% HDI)`, 
                                                  "([0-9].[0-9]{2})e-0([0-9])",
                                                  "$\\1 \times 10^{-\\2}$")
  ) %>%
  dplyr::arrange(Parameter, Population, Garden) %>%
  dplyr::mutate(
    Parameter = dplyr::case_when(
      parameter == "va" ~ "$V_\text{A}$",
      parameter == "vblock" ~ "Block",
      parameter == "vres" ~ "Residual",
      parameter == "h2" ~ "$h^2$"
    )
  ) %>%
  dplyr::select(Population, Garden, Parameter, `Median (95% HDI)`)

# Figure summarizing variance components
df %<>% 
  dplyr::mutate(
    Garden = factor(stringr::str_to_sentence(garden), 
                    levels = c("South", "North")),
    par = factor(dplyr::case_when(
      parameter == "va" ~ "italic(V[A])",
      parameter == "vblock" ~ "Block", 
      parameter == "vres" ~ "Residual", 
      parameter == "h2" ~ "italic(h) ^ 2"
    ), levels = c("Block", "Residual", "italic(V[A])", "italic(h) ^ 2"))
  )

df_summary <- df %>%
  dplyr::group_by(pop, Garden, par) %>%
  tidybayes::point_interval(value, .width = c(0.8, 0.95), .point = median,
                            .interval = tidybayes::hdi) 

ggplot(df, aes(pop, value, fill = Garden)) +
  facet_wrap( ~ par, scales = "free", labeller = label_parsed) +
  geom_violin(trim = FALSE, scale = "width", adjust = 2) +
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
  scale_fill_manual(values = c("grey", "white")) +
  xlab("Population") +
  ylab("Variance or Heritability") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 14),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12)
  )

ggsave("ms/figures/h2-surv.pdf", width = 6.5, height = 6.5, units = "in")

export2ms("vc_table_surv")
