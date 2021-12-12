source("r/header.R")

germination = readr::read_rds("processed-data/germination.rds")
sow_dates = readr::read_rds("r/objects/sow_dates.rds")
census_dates = readr::read_rds("r/objects/census_dates.rds")

file.remove(".gitattributes")
file.create(".gitattributes")
cat("*.rds filter=lfs diff=lfs merge=lfs -text", "\n", file = ".gitattributes")

tidyr::crossing(
  model = c("discrete_lognormal"),
  poph2 = c(FALSE, TRUE)
) %>%
  purrr::pwalk(~{
    stem = glue::glue("stan/{model}_{poph2}.stan", 
                      model = stringr::str_remove(..1, "discrete_"),
                      poph2 = as.numeric(..2))
    write_model(..1, sow_dates, census_dates, germination, stem, nonadd = FALSE,
                poph2 = ..2, interaction = TRUE)
    cat(glue::glue("{stem} filter=lfs diff=lfs merge=lfs -text"), "\n", 
        file = ".gitattributes", append = TRUE)
    cat(glue::glue("{stem} filter=lfs diff=lfs merge=lfs -text\n",
                   stem = stringr::str_replace(stem, "stan", "r/objects") %>%
                     stringr::str_replace(".stan", ".rds")), "\n", 
        file = ".gitattributes", append = TRUE)
  })
