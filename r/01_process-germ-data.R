source("r/header.r")

# Load processed data -----
data = readr::read_csv(
  "raw-data/germination.csv",
  col_types = readr::cols(
    Pop = readr::col_factor(levels = pop_levels()),
    Dam = readr::col_character(),
    Sire = readr::col_character(),
    id = readr::col_character(),
    Cohort = readr::col_character(),
    GermBlock = readr::col_character(),
    Germinate = readr::col_character(),
    SowDate = readr::col_date(format = ""),
    GermDate = readr::col_date(format = ""),
    MainBlock = readr::col_character(),
    InitialSize = readr::col_double(),
    WinterSurv = readr::col_character())
) %>%
  dplyr::mutate(DaysToGerm = as.integer(GermDate - SowDate))

# Identify subcohorts and sow dates -----
sow_dates_north = data %>% 
  dplyr::filter(Cohort == "north") %>% 
  magrittr::use_series(SowDate) %>% 
  unique() %>%
  magrittr::set_names(stringr::str_c("north", 1:length(.)))

sow_dates_south = data %>% 
  dplyr::filter(Cohort == "south") %>% 
  magrittr::use_series(SowDate) %>% 
  unique() %>%
  magrittr::set_names(stringr::str_c("south", 1:length(.)))

sow_dates = tibble::tibble(
  date = c(sow_dates_north, sow_dates_south),
  sub_cohort = c(names(sow_dates_north), names(sow_dates_south))
)

rm(sow_dates_north, sow_dates_south)
readr::write_rds(sow_dates, "r/objects/sow_dates.rds")

# Modify data -----
export2ms("data")

# Remove DaysToGerm < 5 (these are weeds rather than monkeyflowers)
# This also removes NAs
nRemove = length(which(data$DaysToGerm < 5L))
data = dplyr::filter(data, DaysToGerm >= 5L)

# Save modified data -----
readr::write_rds(data, "processed-data/germination.rds")

census_dates = data %>%
  dplyr::group_by(Cohort, GermDate) %>%
  dplyr::summarise() %T>%
  readr::write_rds("r/objects/census_dates.rds")

# Export to ms -----
export2ms(c("nRemove"))
