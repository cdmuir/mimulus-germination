source("r/header.r")

# Load model and calculate additive genetic variance
pars <- c("bPop_germ", "sGeno_germ", "sMat_germ", "sBlock_germ", "sigma")
fit_germ <- readr::read_rds("r/objects/fit.rds")$draws(pars) %>%
  posterior::as_draws_df() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    vpop = var(c(`bPop_germ[1]`, `bPop_germ[2]`, `bPop_germ[3]`, `bPop_germ[4]`, 
                 `bPop_germ[5]`)) * 4/5,
    vg = 4 * sGeno_germ ^ 2,
    vmat = sMat_germ ^ 2,
    vblock = sBlock_germ ^ 2,
    vres = sigma ^ 2,
    vp = vg + vmat + vblock + vres,
    h2 = vg / vp
  ) %>%
  dplyr::ungroup()

# Table summarizing variance components
vc_table_germ <- fit_germ %>%
  dplyr::select(vpop, vg, vmat, vblock, vres, h2) %T>%
  assign(x = "diff_vpop_vg", envir = .GlobalEnv) %>%
  tidyr::pivot_longer(tidyr::everything()) %>%
  dplyr::group_by(name) %>%
  tidybayes::point_interval(.point = median, .interval = tidybayes::qi) %>%
  dplyr::mutate(dplyr::across(value:.upper, signif, digits = 3)) %>%
  dplyr::mutate(
    `Median (95\\% CI)` = glue::glue("${value}~({.lower},~{.upper})$"),
    Parameter = factor(name, levels = c("vpop", "vg", "vmat", "vblock", "vres", "h2"))
  ) %>%
  dplyr::arrange(Parameter) %>%
  dplyr::mutate(
    Parameter = dplyr::case_when(
      name == "vpop" ~ "$V_\\text{pop}$",
      name == "vg" ~ "$V_G$",
      name == "vmat" ~ "$V_M$",
      name == "vblock" ~ "Block",
      name == "vres" ~ "$V_E$",
      name == "h2" ~ "$H^2$"
    )
  ) %>%
  dplyr::select(Parameter, `Median (95\\% CI)`)

# Difference between V_pop and V_G for ms
diff_vpop_vg_germ <- diff_vpop_vg %>%
  dplyr::transmute(diff_vpop_vg = vpop - vg) %>%
  tidybayes::point_interval(.point = median, .interval = tidybayes::qi) 

# Figure summarizing variance components
df <- fit_germ %>%
  dplyr::select(vpop, vg, vmat, vblock, vres, h2) %>%
  tidyr::pivot_longer(tidyr::everything(), names_to = "parameter") %>%
  dplyr::mutate(
    parameter = factor(parameter, levels = c("vpop", "vg", "vmat", "vblock", "vres", "h2")),
    type = dplyr::case_when(
      stringr::str_detect(parameter, "^v[a-z]+$") ~ "Variance components",
      stringr::str_detect(parameter, "^[hH]2$") ~ "Heritability"
  ))

# Variance component panel
gp1 <- ggplot(dplyr::filter(df, type == "Variance components"), 
              aes(parameter, value)) +
  facet_grid(. ~ type) + 
  # geom_violin(trim = FALSE, scale = "width", fill = "grey", adjust = 2) +
  stat_summary(
    fun.data = tidybayes::median_qi, 
    fun.args = list(.width = 0.95),
    geom = "linerange"
  ) +
  stat_summary(
    fun.data = tidybayes::median_qi, 
    fun.args = list(.width = 0.8),
    geom = "linerange",
    size = 1.2
  ) +
  stat_summary(fun = median, geom = "point", size = 3) +
  scale_x_discrete(labels = c(
    expression(italic(V)[pop]), 
    expression(italic(V)[G]), 
    expression(italic(V)[M]), 
    "Block", 
    expression(italic(V)[E])
  )) +
  xlab(element_blank()) +
  ylab("Variance") +
  ggtitle("Germination rate") +
  ylim(0, 0.3) +
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
  # geom_violin(trim = FALSE, scale = "width", fill = "grey", adjust = 2) +
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
  scale_x_discrete(labels = expression(italic(H)^2)) +
  ylim(0, 1) +
  xlab(element_blank()) +
  ylab("Heritability") +
  theme(
    axis.text = element_text(size = 12),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12)
  )

gp <- cowplot::plot_grid(gp1, gp2, nrow = 1, labels = "AUTO", align = "h",
                         rel_widths = c(2/3, 1/3), label_y = 0.95)

ggsave("ms/figures/h2-germ.pdf", plot = gp, width = 6.5, height = 4, 
       units = "in")

export2ms(c("vc_table_germ", "diff_vpop_vg_germ"))
