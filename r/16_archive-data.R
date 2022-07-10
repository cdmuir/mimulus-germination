# Export data file for Dryad
source("r/header.R")

readr::read_rds("processed-data/germination.rds") |>
  dplyr::select(-Germinate) |>
  dplyr::rename(
    population = Pop,
    dam = Dam,
    sire = Sire,
    cohort = Cohort,
    block_germination = GermBlock,
    date_sow = SowDate,
    date_germination = GermDate,
    time_to_germination_days = DaysToGerm,
    block_garden = MainBlock,
    initial_size_mm = InitialSize,
    winter_survival_boolean = WinterSurv
  ) |>
  dplyr::mutate(winter_survival_boolean = (winter_survival_boolean == "YES")) |>
  readr::write_tsv("ms/muir-etal-2022.txt")
