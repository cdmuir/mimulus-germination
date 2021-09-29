# Uncomment if you need to install packages:
# source("r/install-packages.R")

source("r/header.R")

# Uncomment if you need to download large objects from remote repository
# You should only need to do this once
# source("r/00_import-large-objects.R") # slow!

# Emergence model selection
source("r/01_write-germ-models.R")
source("r/02_manipulate-germ-data.R")
source("r/03_extract-and-save.R")
source("r/04_analyze-emergence.R")
source("r/05_check-diagnostics.R")
source("r/06_compare-models.R")

# Comparing additive and nonadditive models
source("r/07_test-nonadditive.R")
source("r/08_prior-vs-posterior.R")
source("r/09_plot-phenotype.R")

# Survival model selection (binomial or beta-binomial)
source("r/10_write-surv-models.R")
source("r/11_manipulate-surv-data.R")
source("r/12_analyze-survival.R")
source("r/13_check-diagnostics.R")
source("r/14_compare-models.R")

# Render manuscript
setwd(path_ms)
rmarkdown::render("ms.Rmd", c("html_document"))
setwd(path)

#knit("ms.Rmd")
#knit("ms.md")
#knit("supplement.Rmd")
