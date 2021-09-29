source("r/header.r")

n_chain = 4

init = readr::read_rds("r/objects/init.rds")

seeds = read.table("r/objects/seeds.txt") %>%
  dplyr::pull(V1)
germ_stan = readr::read_rds("r/objects/germ_stan.rds")
surv_stan = readr::read_rds("r/objects/surv_stan.rds")
data = c(germ_stan, surv_stan)

# Remove unneeded elements from data that have NAs or are character
vars_to_drop = c(
  tidyr::crossing(
    x = c("Dam", "Sire", "SireDam"),
    pop = pop_levels(),
    y = c("_germ", "_surv")
  ) %>%
    dplyr::transmute(z = paste0(x, pop, y)) %>%
    dplyr::pull(z),
  purrr::map_lgl(data, ~{!is.numeric(.x)}) %>%
    which() %>%
    names()
)
data = data[setdiff(names(data), vars_to_drop)]

data$ntd = rep(4, data$nPop_germ)
# check that NAs removed
# which(purrr::map_lgl(data, function(.x) any(is.na(.x))))

# Fit models -----
models = c(
  "lognormal_0_0_0",
  "lognormal_0_1_0",
  "lognormal_0_0_1",
  "lognormal_0_1_1"
)

# for (i in seq_along(models)) {
i = 3
  mod = cmdstan_model(glue::glue("stan/{stem}.stan", stem = models[i]))
  
  fit = mod$sample(
    data = data,
    seed = seeds[i],
    chains = n_chain, 
    parallel_chains = n_chain,
    refresh = 4e2,
    iter_warmup = 4e3,
    iter_sampling = 4e3,
    thin = 4e0,
    max_treedepth = 10L
  )

  fit$save_object(glue::glue("r/objects/{stem}.rds", stem = models[i]))
  
}
