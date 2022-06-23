source("r/header.R")

world = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
states = sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

card_pops = readr::read_csv("raw-data/card_populations.csv") |>
  dplyr::select(longitude = Lon, latitude = Lat) |>
  dplyr::filter(longitude < -113)

focal_pops = tibble::tribble(
    ~Name, ~pop, ~latitude, ~longitude, ~`Elevation (mas)`, ~garden,
    "Sweetwater River", "CUR", 32.900, -116.585, 1180, "bold",
    "West Fork Mojave River", "WFM", 34.284, -117.378, 1120, "plain",
    "North Fork Middle Tule River", "NMT", 36.201, -118.759, 926, "plain",
    "Little Jamison Creek", "LIJ", 39.743, -120.704, 1603, "bold",
    "Rock Creek", "RCK", 43.374, -122.957, 326, "plain"
  ) |>
  dplyr::mutate(pop = factor(pop, levels = rev(pop_levels())))

gp = ggplot(data = world) +
  geom_sf() +
  geom_sf(data = states, fill = NA) + 
  geom_point(data = card_pops, aes(x = longitude, y = latitude), 
                      size = 4, fill = "darkgrey") +
  geom_point(data = focal_pops, aes(x = longitude, y = latitude, color = pop), 
                      size = 4, show.legend = FALSE ) +
  geom_label(
    data = focal_pops, 
    aes(x = longitude + c(2, 2.5, 2.5, 2.25, 1.5), y = latitude, 
        fontface = garden, fill = garden, label = pop_levels() %>% 
          sapply(get_labels) %>% sapply(place_line_break)),
    hjust = 0.5, show.legend = FALSE
  ) +
  scale_color_manual(values = rev(palette())) +
  scale_fill_manual(values = c("white", "lightgrey")) +
  coord_sf(xlim = c(-125, -113), ylim = c(32, 44), expand = FALSE) +
  ggtitle("Source populations within species' range")

ggsave("ms/figures/range-map.pdf", plot = gp, width = 5.25, height = 7, units = "in")
