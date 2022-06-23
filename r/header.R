rm(list = ls())
graphics.off()

# Libraries
library(cmdstanr)
library(cowplot)
library(ggplot2)
library(magrittr)

source("r/functions.R")

seeds = readr::read_lines("raw-data/seeds.txt")

palette(colorRampPalette(c("tomato", "steelblue"), alpha = TRUE)(5))

theme_set(theme_cowplot())