source("r/header.r")

# Load modified data -----
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
)

message("just remove NAs right now for testing model. need to recheck.")
data %<>% dplyr::filter(!is.na(WinterSurv))

# Convert survival to integer
data$WinterSurv %<>%
  stringr::str_replace("NO", "0") %>%
  stringr::str_replace("YES", "1") %>%
  as.integer()
  
# Manipulate data for Stan -----
surv_stan <- data %>%
  dplyr::select(Pop_surv = Pop, Dam_surv = Dam, Sire_surv = Sire, 
                Cohort_surv = Cohort, Block_surv = MainBlock, WinterSurv) %>%
  dplyr::filter(complete.cases(.)) %>%
  as.list()

surv_stan$Cohort_surv %<>%
  stringr::str_replace("south", "0") %>%
  stringr::str_replace("north", "1") %>%
  as.integer()

surv_stan$Block_surv %<>% as.factor()
surv_stan$nBlock_surv <- nlevels(surv_stan$Block_surv)
surv_stan$Block_surv %<>% as.integer()
surv_stan$nPop_surv <- nlevels(surv_stan$Pop_surv)
  
for (pop in levels(surv_stan$Pop_surv)) {
    
  nObs <- length(surv_stan$Pop_surv)
  x <- surv_stan$Pop_surv == pop
    
  surv_stan[[stringr::str_c("Sire", pop, "_surv")]] <- 
    surv_stan[[stringr::str_c("Dam", pop, "_surv")]] <-
    surv_stan[[stringr::str_c("SireDam", pop, "_surv")]] <- rep(NA, nObs)
    
  surv_stan[[stringr::str_c("Sire", pop)]][x] <- 
    subset(surv_stan$Sire_surv, x)
  surv_stan[[stringr::str_c("Dam", pop)]][x] <- 
    subset(surv_stan$Dam_surv, x)
  surv_stan[[stringr::str_c("SireDam", pop)]][x] <- 
    stringr::str_c(surv_stan[[stringr::str_c("Sire", pop)]][x], 
                   surv_stan[[stringr::str_c("Dam", pop)]][x],
                   sep = "_")
  
  surv_stan[[stringr::str_c("Sire", pop)]] %<>% as.factor() %>% as.integer()
  surv_stan[[stringr::str_c("Dam", pop)]] %<>% as.factor() %>% as.integer()
  surv_stan[[stringr::str_c("SireDam", pop)]] %<>% as.factor() %>% as.integer()
    
  surv_stan[[stringr::str_c("nSire", pop, "_surv")]] <- 
    surv_stan[[stringr::str_c("Sire", pop)]] %>% 
    na.omit() %>% unique() %>% length()
  surv_stan[[stringr::str_c("nDam", pop, "_surv")]] <- 
    surv_stan[[stringr::str_c("Dam", pop)]] %>%
    na.omit() %>% unique() %>% length()
  surv_stan[[stringr::str_c("nSireDam", pop, "_surv")]] <- 
    surv_stan[[stringr::str_c("SireDam", pop)]] %>%
    na.omit() %>% unique() %>% length()
    
  stopifnot(surv_stan[[stringr::str_c("nSire", pop, "_surv")]] == surv_stan[[stringr::str_c("nDam", pop, "_surv")]])
 
  surv_stan[[stringr::str_c("n", pop)]] <- length(which(x))
  surv_stan[[stringr::str_c("Cohort", pop)]] <- surv_stan$Cohort_surv[x]
  surv_stan[[stringr::str_c("Block", pop)]] <- surv_stan$Block_surv[x]
  surv_stan[[stringr::str_c("Sire", pop)]] <- 
    surv_stan[[stringr::str_c("Sire", pop)]][x]
  surv_stan[[stringr::str_c("Dam", pop)]] <- 
    surv_stan[[stringr::str_c("Dam", pop)]][x]
  surv_stan[[stringr::str_c("SireDam", pop)]] <- 
    surv_stan[[stringr::str_c("SireDam", pop)]][x]
  surv_stan[[stringr::str_c("WinterSurv", pop)]] <- 
    surv_stan$WinterSurv[x]
  
}
  
# Save manipulated dataset -----
readr::write_rds(surv_stan, "r/objects/surv_stan.rds")
