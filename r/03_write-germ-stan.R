source("r/header.r")

# Load modified data -----
data = readr::read_rds("processed-data/germination.rds")

# Load sow dates
sow_dates = readr::read_rds("r/objects/sow_dates.rds") %>% 
  dplyr::mutate(SowDate = date, SubCohort = sub_cohort) %>% 
  dplyr::select(-date) 

data %<>% dplyr::full_join(sow_dates, by = "SowDate")

# Manipulate data for Stan -----
germ_stan = dplyr::select(data, Pop_germ = Pop, Dam_germ = Dam, Sire_germ = Sire, 
                          Cohort_germ = Cohort, SubCohort, 
                          Block_germ = GermBlock, SowDate, GermDate)
germ_stan = germ_stan[complete.cases(germ_stan), ]
germ_stan %<>% 
  dplyr::mutate(DaysToGerm = as.integer(GermDate - SowDate)) %>% 
  as.list()
germ_stan$Cohort_germ = ifelse(germ_stan$Cohort_germ == "north", 1, 0)
germ_stan$Block_germ %<>% as.factor()
germ_stan$nBlock_germ = nlevels(germ_stan$Block_germ)
germ_stan$Block_germ %<>% as.integer()
germ_stan$Pop_germ %<>% as.factor()
germ_stan$nPop_germ = nlevels(germ_stan$Pop_germ)

for (pop in levels(germ_stan$Pop_germ)) {
  nObs_germ = length(germ_stan$Pop_germ)
  x = germ_stan$Pop_germ == pop
    
  # need to change to length == nObs, all from other pops are NA
  germ_stan[[stringr::str_c("Sire", pop, "_germ")]] = 
    germ_stan[[stringr::str_c("Dam", pop, "_germ")]] =
    germ_stan[[stringr::str_c("SireDam", pop, "_germ")]] = rep(NA, nObs_germ)
  
  germ_stan[[stringr::str_c("Sire", pop, "_germ")]][x] = subset(germ_stan$Sire_germ, x)
  germ_stan[[stringr::str_c("Dam", pop, "_germ")]][x] = subset(germ_stan$Dam_germ, x)
  germ_stan[[stringr::str_c("SireDam", pop, "_germ")]][x] = 
    stringr::str_c(germ_stan[[stringr::str_c("Sire", pop)]][x], 
                   germ_stan[[stringr::str_c("Dam", pop)]][x],
                   "germ", sep = "_")
  
  germ_stan[[stringr::str_c("Sire", pop, "_germ")]] %<>% as.factor() %>% as.integer()
  germ_stan[[stringr::str_c("Dam", pop, "_germ")]] %<>% as.factor() %>% as.integer()
  germ_stan[[stringr::str_c("SireDam", pop, "_germ")]] %<>% as.factor() %>% as.integer()
  
  germ_stan[[stringr::str_c("nSire", pop, "_germ")]] = 
    germ_stan[[stringr::str_c("Sire", pop, "_germ")]] %>%
    na.omit() %>% unique() %>% length()
  germ_stan[[stringr::str_c("nDam", pop, "_germ")]] = 
    germ_stan[[stringr::str_c("Dam", pop, "_germ")]] %>%
    na.omit() %>% unique() %>% length()
  germ_stan[[stringr::str_c("nSireDam", pop, "_germ")]] = 
    germ_stan[[stringr::str_c("SireDam", pop, "_germ")]] %>%
    na.omit() %>% unique() %>% length()
  
  stopifnot(germ_stan[[stringr::str_c("nSire", pop, "_germ")]] == 
              germ_stan[[stringr::str_c("nDam", pop, "_germ")]])
  
}

for (sub_cohort in sow_dates$SubCohort) {
  
  for (pop in levels(germ_stan$Pop_germ)) {
    
    cohort = ifelse(stringr::str_detect(sub_cohort, "north"), "north", "south")
    suffix = stringr::str_c("_", sub_cohort, "_", pop)
    
    x = which(germ_stan$SubCohort == sub_cohort & germ_stan$Pop_germ == pop)
    
    germ_stan[[stringr::str_c("n", suffix)]] = length(x)
    germ_stan[[stringr::str_c("Cohort", suffix)]] = germ_stan$Cohort_germ[x]
    germ_stan[[stringr::str_c("Block", suffix)]] = germ_stan$Block_germ[x]
    germ_stan[[stringr::str_c("Sire", suffix)]] = 
      germ_stan[[stringr::str_c("Sire", pop, "_germ")]][x]
    germ_stan[[stringr::str_c("Dam", suffix)]] = 
      germ_stan[[stringr::str_c("Dam", pop, "_germ")]][x]
    germ_stan[[stringr::str_c("SireDam", suffix)]] = 
      germ_stan[[stringr::str_c("SireDam", pop, "_germ")]][x]
    germ_stan[[stringr::str_c("DaysToGerm", suffix)]] = 
      germ_stan$DaysToGerm[x]
    
  }
}

# Save manipulated dataset -----
readr::write_rds(germ_stan, "r/objects/germ_stan.rds")
