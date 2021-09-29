# A good simple tutorial about Make can be found at http://kbroman.org/minimal_make/ 
R_OPTS=--no-save --no-restore --no-init-file --no-site-file
# sections copied from other project are commented out. modify and uncomment when readonly
# all: data model paper
# data: objects/pair_div.rds 
# model: objects/fit.rds
# paper: ms/ms.pdf figures/concepts.pdf figures/h1-raw.pdf figures/h1.pdf figures/h2.pdf

# 01_process-germ-data.R ----
r/objects/sow_dates.rds: raw-data/germination.csv r/01_process-germ-data.R
	Rscript -e 'source("r/01_process-germ-data.R")'

processed-data/germination.rds: raw-data/germination.csv r/01_process-germ-data.R
	Rscript -e 'source("r/01_process-germ-data.R")'

r/objects/census_dates.rds: raw-data/germination.csv r/01_process-germ-data.R
	Rscript -e 'source("r/01_process-germ-data.R")'

# 02_write-germ-stan.R ----
r/objects/germ_stan.rds: r/objects/sow_dates.rds processed-data/germination.rds r/02_write-germ-stan.R
	Rscript -e 'source("r/02_write-germ-stan.R")'

# 03_write-surv-stan.R ----
r/objects/surv_stan.rds: raw-data/germination.csv r/03_write-surv-stan.R
	Rscript -e 'source("r/03_write-surv-stan.R")'

# 04_write-models.R ----
.gitattributes: r/04_write-models.R
	Rscript -e 'source("r/04_write-models.R")'
stan/lognormal_0_0_0.stan: processed-data/germination.rds r/objects/sow_dates.rds r/objects/census_dates.rds r/04_write-models.R
	Rscript -e 'source("r/04_write-models.R")'
stan/lognormal_0_0_1.stan: processed-data/germination.rds r/objects/sow_dates.rds r/objects/census_dates.rds r/04_write-models.R
	Rscript -e 'source("r/04_write-models.R")'
stan/lognormal_0_1_0.stan: processed-data/germination.rds r/objects/sow_dates.rds r/objects/census_dates.rds r/04_write-models.R
	Rscript -e 'source("r/04_write-models.R")'
stan/lognormal_0_1_1.stan: processed-data/germination.rds r/objects/sow_dates.rds r/objects/census_dates.rds r/04_write-models.R
	Rscript -e 'source("r/04_write-models.R")'
stan/lognormal_1_0_0.stan: processed-data/germination.rds r/objects/sow_dates.rds r/objects/census_dates.rds r/04_write-models.R
	Rscript -e 'source("r/04_write-models.R")'
stan/lognormal_1_0_1.stan: processed-data/germination.rds r/objects/sow_dates.rds r/objects/census_dates.rds r/04_write-models.R
	Rscript -e 'source("r/04_write-models.R")'
stan/lognormal_1_1_0.stan: processed-data/germination.rds r/objects/sow_dates.rds r/objects/census_dates.rds r/04_write-models.R
	Rscript -e 'source("r/04_write-models.R")'
stan/lognormal_1_1_1.stan: processed-data/germination.rds r/objects/sow_dates.rds r/objects/census_dates.rds r/04_write-models.R
	Rscript -e 'source("r/04_write-models.R")'

# 05_fit-models.R ----
r/objects/lognormal_0_0_0.rds: stan/lognormal_0_0_0.stan r/05_fit-models.R
	Rscript -e 'source("r/05_fit-models.R")'
r/objects/lognormal_0_0_1.rds: stan/lognormal_0_0_1.stan r/05_fit-models.R
	Rscript -e 'source("r/05_fit-models.R")'
r/objects/lognormal_0_1_0.rds: stan/lognormal_0_1_0.stan r/05_fit-models.R
	Rscript -e 'source("r/05_fit-models.R")'
r/objects/lognormal_0_1_1.rds: stan/lognormal_0_1_1.stan r/05_fit-models.R
	Rscript -e 'source("r/05_fit-models.R")'
r/objects/lognormal_1_0_0.rds: stan/lognormal_1_0_0.stan r/05_fit-models.R
	Rscript -e 'source("r/05_fit-models.R")'
r/objects/lognormal_1_0_1.rds: stan/lognormal_1_0_1.stan r/05_fit-models.R
	Rscript -e 'source("r/05_fit-models.R")'
r/objects/lognormal_1_1_0.rds: stan/lognormal_1_1_0.stan r/05_fit-models.R
	Rscript -e 'source("r/05_fit-models.R")'
r/objects/lognormal_1_1_1.rds: stan/lognormal_1_1_1.stan r/05_fit-models.R
	Rscript -e 'source("r/05_fit-models.R")'

# copied from other project
# ms/ms.pdf: ms/ms.Rmd ms/stomata-independence.bib figures/concepts.pdf figures/h1-raw.pdf figures/h2-raw.pdf figures/h1.pdf figures/h2.pdf objects/modeloutput.rds
#	Rscript -e 'rmarkdown::render("ms/ms.Rmd", output_format = "bookdown::pdf_document2", output_file = "ms.pdf")'

figures/concepts.pdf: r/16_plot-concepts.R
	Rscript -e 'source("r/16_plot-concepts.R")'

figures/h1-raw.pdf: r/08_plot-h1-raw.R
	Rscript -e 'source("r/08_plot-h1-raw.R")'

figures/h2-raw.pdf: r/09_plot-h2-raw.R
	Rscript -e 'source("r/09_plot-h2-raw.R")'

figures/h1.pdf: r/11_prepare-plotting.R r/12_plot-h1-pairs.R r/14_plot-h1.R objects/fit.rds
	Rscript -e 'source("r/11_prepare-plotting.R")'
	Rscript -e 'source("r/12_plot-h1-pairs.R")'
	Rscript -e 'source("r/14_plot-h1.R")'

figures/h2.pdf: r/11_prepare-plotting.R r/13_plot-h2-pairs.R r/15_plot-h2.R objects/fit.rds
	Rscript -e 'source("r/11_prepare-plotting.R")'
	Rscript -e 'source("r/13_plot-h2-pairs.R")'
	Rscript -e 'source("r/15_plot-h2.R")'

objects/modeloutput.rds: objects/fit.rds r/17_summarize-fit.R
	Rscript -e 'source("r/17_summarize-fit.R")'

objects/fit.rds: r/10_fit-pairs.R processed-data/trimmed-data.rds processed-data/trimmed-phylogeny.rds objects/pair_div.rds
	Rscript -e 'source("r/10_fit-pairs.R")'

processed-data/trimmed-data.rds: r/06_connect-data.R processed-data/resolved_names.rds processed-data/phy2.rds
	Rscript -e 'source("r/06_connect-data.R")'
	
processed-data/trimmed-phylogeny.rds: r/06_connect-data.R processed-data/resolved_names.rds processed-data/phy2.rds
	Rscript -e 'source("r/06_connect-data.R")'

objects/pair_div.rds: r/07_make-pairs.R	processed-data/trimmed-data.rds processed-data/trimmed-phylogeny.rds
	Rscript -e 'source("r/07_make-pairs.R")'

processed-data/resolved_names.rds: r/02_resolve-names.R processed-data/taxize_output.rds
	Rscript -e 'source("r/02_resolve-names.R")'

processed-data/phy2.rds: r/05_make-subtrees.R processed-data/phy1.rds objects/polytomies.rds
	Rscript -e 'source("r/05_make-subtrees.R")'
	
processed-data/taxize_output.rds: r/01_taxize-data.R
	Rscript -e 'source("r/01_taxize-data.R")'

processed-data/phy1.rds: r/04_prepare-subtrees.R processed-data/phy.rds
	Rscript -e 'source("r/04_prepare-subtrees.R")'
	source raxml/run-raxml.sh

objects/polytomies.rds: r/04_prepare-subtrees.R processed-data/phy.rds
	Rscript -e 'source("r/04_prepare-subtrees.R")'

processed-data/phy.rds: r/03_make-megatree.R processed-data/resolved_names.rds
	Rscript -e 'source("r/03_make-megatree.R")'

clean: 
	\rm -f *~ *.Rout */*~ */*.Rout .RData Rplots.pdf
	
cleanall: 
	\rm -f *.aux *.bbl *.blg *.log *.pdf *~ *.Rout */*~ */*.Rout figures/*.png figures/*.pdf ms/ms.pdf ms/ms.tex objects/*.rds processed-data/*.rds */*.aux */*.log 

