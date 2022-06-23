source("r/header.R")

# Data accessed: 2021-11-29 from http://www.climatewna.com/default.aspx
# Variable descriptions: http://www.climatewna.com/help/ClimateNA/Help.htm

# Seasons:
# Winter (_wt): Dec. (prev. yr) - Feb for annual, Jan, Feb, Dec for normals
# Spring (_sp): March, April and May
# Summer (_sm): June, July and August
# Autumn (_at): September, October and November

climate_data = readr::read_csv("raw-data/climate_data.csv") |>
  dplyr::select(pop, garden, period, Tave_at, Tave_wt, Tave_sp) |>
  tidyr::pivot_longer(Tave_at:Tave_sp) |>
  dplyr::filter(!(
    (period == "Year_2015.ann" & name == "Tave_wt") |
      (period == "Year_2015.ann" & name == "Tave_sp") |
      (period == "Year_2016.ann" & name == "Tave_at")
  )) |>
  dplyr::mutate(
    pop = factor(pop, levels = pop_levels()),
    Garden = factor(stringr::str_to_sentence(garden), 
                    levels = c("South", "North")),
    season = dplyr::case_when(
      name == "Tave_at" ~ 1,
      name == "Tave_wt" ~ 2,
      name == "Tave_sp" ~ 3
    ))

gp <- ggplot() +
  geom_line(
    data = dplyr::filter(climate_data, !is.na(Garden)),
    mapping = aes(season, value, linetype = Garden),
    size = 1.5, lineend = "round"
  ) +
  geom_point(
    data = dplyr::filter(climate_data, !is.na(pop)), 
    mapping = aes(season, value, color = pop), size = 3
  ) +
  scale_color_manual(
    values = palette(),
    labels = pop_levels() %>% sapply(get_labels) %>% sapply(place_line_break),
    name = "Source\nPopulation"
  ) +
  scale_x_continuous(breaks = 1:3, labels = c("autumn", "winter", "spring"),
                     limits = c(0.75, 3.25)) +
  xlab("Season") +
  ylab(expression(Average~temperature~(degree~C))) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.key.width = unit(1, "cm"),
    legend.position = "right",
    legend.text = element_text(size = 12, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(size = 14, hjust = 0.5)
  )

ggsave("ms/figures/climate.pdf", plot = gp, width = 6.5, height = 4, 
       units = "in")
