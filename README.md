# mimulus-germination
Are phenotypic clines adaptive?

This project was developed by [Chris Muir](https://cdmuir.netlify.app) and [Amy Angert](https://www.botany.ubc.ca/people/amy-angert/). [Courtney van den Elzen](https://www.colorado.edu/lab/emery/courtney-van-den-elzen) helped with data collection.

More information about the study is available in a preprint which you can find on [biorxiv](https://doi.org/10.1101/######) or on [github](https://github.com/cdmuir/mimulus-germination/blob/master/ms/ms.pdf).

## Instructions to reproduce manuscript using RStudio

While you can run the source code for this manuscript through any R interface, I strongly recommend using the [RStudio](https://www.rstudio.com/) project file included in this repository.

1. Download or clone this repository to your machine. To clone from the Terminal:

```
git clone git@github.com:cdmuir/mimulus-germination.git
```

2. Open `mimulus-germination.Rproj` in [RStudio](https://www.rstudio.com/)

``` {r}
source("r/install-packages.R")
```

3. Install R packages if necessary. Running "r/install-packages.R" will do this for you.

If you don't want to rerun analyses (slow) or refit Stan models (*very* slow), you are ready to generate the manuscript on your machine. If you want to rerun analyses or refit models, see instructions below.

## Generating manuscript

Simply open `ms/ms.Rnw` and compile using RStudio.

## Running analyses (slow)

1. Open `r/run-all.R` to begin running scripts. 

2. Import large objects (mostly posterior samples from Stan models). From the R Console, do:

``` {r}
source("r/00_import-large-objects.R") 
```
