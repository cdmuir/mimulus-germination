rm(list = ls())
graphics.off()

# Libraries
library(cmdstanr)
library(cowplot)
library(ggplot2)
library(magrittr)

if (file.exists("r/private.R")) source("r/private.R")

source("r/functions.R")

seeds = readr::read_lines("r/objects/seeds.txt")

palette(colorRampPalette(c("tomato", "steelblue"), alpha = TRUE)(5))
