source("r/header.r")

# SCRATCH TO HELP calculating H2 for lognormal trait

N = 1e6; set.seed(1234); XX = rlnorm(N, rnorm(N, 0, sqrt(4)), sqrt(10))
mean(XX)
var(XX)

mean(log(XX)) # should be 0
var(log(XX)) # should be 14

# Load model and calculate additive genetic variance
pars <- c("bPop_germ", "sGeno_germ", "sBlock_germ", "sigma")
fit_germ <- readr::read_rds("r/objects/fit.rds")$draws(pars) %>%
  posterior::as_draws_df() %>%
  dplyr::mutate(
    va = 4 * sGeno_germ ^ 2,
    vblock = sBlock_germ ^ 2,
    vres = sigma ^ 2,
    vp = va + vblock + vres,
    h2 = va / vp
  )

# Table summarizing variance components
vc_table_germ <- fit_germ %>%
  dplyr::select(va, vblock, vres, h2) %>%
  tidyr::pivot_longer(tidyr::everything()) %>%
  dplyr::group_by(name) %>%
  tidybayes::point_interval(.point = median, .interval = tidybayes::qi) %>%
  dplyr::mutate(dplyr::across(value:.upper, signif, digits = 3)) %>%
  dplyr::mutate(
    `Median (95% CI)` = glue::glue("{.point} ({.lower}--{.upper})"),
    Parameter = factor(name, levels = c("va", "vblock", "vres", "h2"))
  ) %>%
  dplyr::arrange(Parameter) %>%
  dplyr::mutate(
    Parameter = dplyr::case_when(
      name == "va" ~ "$V_\text{A}$",
      name == "vblock" ~ "Block",
      name == "vres" ~ "Residual",
      name == "h2" ~ "$h^2$"
    )
  ) %>%
  dplyr::select(Parameter, `Median (95% CI)`)

# Figure summarizing variance components
df <- fit_germ %>%
  dplyr::select(va, vblock, vres, h2) %>%
  tidyr::pivot_longer(tidyr::everything(), names_to = "parameter") %>%
  dplyr::mutate(
    parameter = factor(parameter, levels = c("va", "vblock", "vres", "h2")),
    type = dplyr::case_when(
      stringr::str_detect(parameter, "^v[a-z]+$") ~ "Variance components",
      stringr::str_detect(parameter, "^[hH]2$") ~ "Heritability"
  ))

# Variance component panel
gp1 <- ggplot(dplyr::filter(df, type == "Variance components"), 
              aes(parameter, value)) +
  facet_grid(. ~ type) + 
  geom_violin(trim = FALSE, scale = "width", fill = "grey", adjust = 2) +
  stat_summary(
    fun.data = tidybayes::median_hdi, 
    fun.args = list(.width = 0.95),
    geom = "linerange"
  ) +
  stat_summary(
    fun.data = tidybayes::median_hdi, 
    fun.args = list(.width = 0.8),
    geom = "linerange",
    size = 1.2
  ) +
  stat_summary(fun = median, geom = "point", size = 3) +
  scale_x_discrete(labels = c(
    expression(italic(V[A])), 
    "Block", 
    "Residual"
  )) +
  xlab(element_blank()) +
  ylab("Variance") +
  theme_cowplot() +
  theme(
    axis.text = element_text(size = 12),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12)
  )
  
# Heritability panel
gp2 <- ggplot(dplyr::filter(df, type == "Heritability"), 
              aes(parameter, value)) +
  facet_grid(. ~ type) + 
  geom_violin(trim = FALSE, scale = "width", fill = "grey", adjust = 2) +
  stat_summary(
    fun.data = tidybayes::median_hdi, 
    fun.args = list(.width = 0.95),
    geom = "linerange"
  ) +
  stat_summary(
    fun.data = tidybayes::median_hdi, 
    fun.args = list(.width = 0.8),
    geom = "linerange",
    size = 1.2
  ) +
  stat_summary(fun = median, geom = "point", size = 3) +
  scale_x_discrete(labels = "Narrow\nsense") +
  ylim(0, 1) +
  xlab(element_blank()) +
  ylab("Heritability") +
  theme_cowplot() +
  theme(
    axis.text = element_text(size = 12),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12)
  )

cowplot::plot_grid(gp1, gp2, nrow = 1, labels = "AUTO", align = "h")

ggsave("ms/figures/h2-germ.pdf", width = 6.5, height = 4, units = "in")

export2ms("vc_table_germ")
