# mimulus-germination

This repository contains source code associated with the manuscript:

[Selection on early survival does not explain germination rate clines in *Mimulus cardinalis*](https://doi.org/10.1002/ajb2.XXXXX). Accepted at *American Journal of Botany*.

This project was developed by [Chris Muir](https://cdmuir.netlify.app) and [Amy Angert](https://www.botany.ubc.ca/people/amy-angert/). [Courtney Van Den Elzen](https://www.colorado.edu/lab/emery/courtney-van-den-elzen) helped with data collection.

## Contents

This repository has the following file folders:

- `ms`: manuscript input (e.g. `ms.Rmd` and `mimulus-germination.bib`) and output (e.g. `ms.pdf` and `si.pdf`) files
- `ms/figures`: figures generated from *R* code
- `ms/objects`: saved objects generated from *R* code
- `processed-data`: processed data generated from *R* code
- `r`: *R* scripts for all data processing and analysis
- `raw-data`: raw data files
- `stan`: [*Stan*](https://mc-stan.org) code for statistical models
- `template`: manuscript style files

## Prerequisites:

To run code and render manuscript:

- [*R*](https://cran.r-project.org/) version >4.1.0 and [*RStudio*](https://www.rstudio.com/)
- [LaTeX](https://www.latex-project.org/): you can install the full version or try [**tinytex**](https://yihui.org/tinytex/)
- [GNU Make](https://www.gnu.org/software/make/): In terminal, you can just type `make paper` to render the manuscript. You can also use it to re-run all scripts.

Before running scripts, you'll need to install the following *R* packages:

```
source("r/install-packages.R")
```

To fit *Stan* models, set up [**cmdstanr**](https://mc-stan.org/cmdstanr/).

## Downloading data and code 

1. Download or clone this repository to your machine.

```
git clone git@github.com:cdmuir/mimulus-germination.git
```

2. Open `mimulus-germination.Rproj` in [RStudio](https://www.rstudio.com/)

## Rendering manuscript

### Software requirements

At minimum, you will need [R](https://cran.r-project.org/) installed on your machine. Install additional packages by running `r/install-packages.R`.

### Rendering manuscript with pre-saved outout

Knit `ms/ms.Rmd`, `ms/tables-andcaptions.Rmd`, and `ms/si.Rmd` using [RStudio](https://www.rstudio.com/).

You can also run the following code from the R console:

```{r}
# Main paper
rmarkdown::render(
  input = "ms/ms.Rmd",
  output_dir = "ms"
)

# Supporting information
rmarkdown::render(
  input = "ms/si.Rmd",
  output_dir = "ms"
)
```

or use `make`

```
make paper
```

### Generating all results

You can re-run all analyses, figures, etc. using [GNU make](https://www.gnu.org/software/make/). Type `make --version` in a terminal/shell to see if it is already installed.

```
# Clear out previously saved output
make cleanall
# This will take a long time to run
make
```
