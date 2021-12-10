source("r/header.R")

germination = readr::read_rds("processed-data/germination.rds")
sow_dates = readr::read_rds("r/objects/sow_dates.rds")
census_dates = readr::read_rds("r/objects/census_dates.rds")

file.remove(".gitattributes")
file.create(".gitattributes")
cat("*.rds filter=lfs diff=lfs merge=lfs -text", "\n", file = ".gitattributes")

tidyr::crossing(
  model = c("discrete_lognormal"),
  nonadd = c(FALSE, TRUE),
  poph2 = c(FALSE, TRUE),
  interaction = c(FALSE, TRUE)
) %>%
  purrr::pwalk(~{
    stem = glue::glue("stan/{model}_{nonadd}_{poph2}_{interaction}.stan", 
                      model = stringr::str_remove(..1, "discrete_"),
                      nonadd = as.numeric(..2), 
                      poph2 = as.numeric(..3), 
                      interaction = as.numeric(..4))
    write_model(..1, sow_dates, census_dates, germination, stem,
                nonadd = ..2, poph2 = ..3, interaction = ..4)
    cat(glue::glue("{stem} filter=lfs diff=lfs merge=lfs -text"), "\n", 
        file = ".gitattributes", append = TRUE)
    cat(glue::glue("{stem} filter=lfs diff=lfs merge=lfs -text\n",
                   stem = stringr::str_replace(stem, "stan", "r/objects") %>%
                     stringr::str_replace(".stan", ".rds")), "\n", 
        file = ".gitattributes", append = TRUE)
  })
