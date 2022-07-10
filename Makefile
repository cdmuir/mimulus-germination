# A good simple tutorial about Make can be found at http://kbroman.org/minimal_make/ 
R_OPTS=--no-save --no-restore --no-init-file --no-site-file

all: data model paper
data: processed-data/germination.rds 
model: r/objects/fit.rds
paper: ms/ms.pdf ms/si.pdf ms/export/data.rds ms/export/df_ind_lm.rds ms/export/df_pop_lm.rds ms/export/diff_vpop_vg_germ.rds ms/export/nRemove.rds ms/export/vc_table_germ.rds ms/export/vc_table_surv.rds ms/export/diff_vpop_vg_germ.rds ms/figures/climate.pdf ms/figures/h2-germ.pdf ms/figures/h2-surv.pdf ms/figures/mean-traits.pdf ms/figures/pp_check_germ.pdf ms/figures/selection.pdf

# 01_process-germ-data.R ----
r/objects/sow_dates.rds: raw-data/germination.csv r/functions.R r/header.R r/01_process-germ-data.R
	Rscript -e 'source("r/01_process-germ-data.R")'
processed-data/germination.rds: raw-data/germination.csv r/functions.R r/header.R r/01_process-germ-data.R
	Rscript -e 'source("r/01_process-germ-data.R")'
r/objects/census_dates.rds: raw-data/germination.csv r/functions.R r/header.R r/01_process-germ-data.R
	Rscript -e 'source("r/01_process-germ-data.R")'
ms/export/data.rds: raw-data/germination.csv r/functions.R r/header.R r/01_process-germ-data.R
	Rscript -e 'source("r/01_process-germ-data.R")'
ms/export/nRemove.rds: raw-data/germination.csv r/functions.R r/header.R r/01_process-germ-data.R
	Rscript -e 'source("r/01_process-germ-data.R")'

# 02_write-models.R ----
.gitattributes: r/02_write-models.R
	Rscript -e 'source("r/02_write-models.R")'
stan/lognormal_0_1.stan: processed-data/germination.rds r/objects/sow_dates.rds r/objects/census_dates.rds r/functions.R r/header.R r/02_write-models.R
	Rscript -e 'source("r/02_write-models.R")'
stan/lognormal_1_1.stan: processed-data/germination.rds r/objects/sow_dates.rds r/objects/census_dates.rds r/functions.R r/header.R r/02_write-models.R
	Rscript -e 'source("r/02_write-models.R")'

# 03_write-germ-stan.R ----
r/objects/germ_stan.rds: r/objects/sow_dates.rds processed-data/germination.rds r/functions.R r/header.R r/03_write-germ-stan.R
	Rscript -e 'source("r/03_write-germ-stan.R")'

# 04_write-surv-stan.R ----
r/objects/surv_stan.rds: raw-data/germination.csv r/functions.R r/header.R r/04_write-surv-stan.R
	Rscript -e 'source("r/04_write-surv-stan.R")'

# 05_fit-models.R ----
r/objects/lognormal_0_1.rds: stan/lognormal_0_1.stan raw-data/seeds.txt r/objects/germ_stan.rds r/objects/surv_stan.rds r/functions.R r/header.R r/05_fit-models.R
	Rscript -e 'source("r/05_fit-models.R")'
r/objects/lognormal_1_1.rds: stan/lognormal_1_1.stan raw-data/seeds.txt r/objects/germ_stan.rds r/objects/surv_stan.rds r/functions.R r/header.R r/05_fit-models.R
	Rscript -e 'source("r/05_fit-models.R")'

# 06_compare-models.R ----
r/objects/fit.rds: r/objects/lognormal_0_1.rds r/objects/lognormal_1_1.rds r/functions.R r/header.R r/06_compare-models.R
	Rscript -e 'source("r/06_compare-models.R")'

# 08_plot-germ-pp.R ----
ms/figures/pp_check_germ.pdf: r/objects/fit.rds r/objects/sow_dates.rds r/objects/germ_stan.rds r/functions.R r/header.R r/08_plot-germ-pp.R
	Rscript -e 'source("r/08_plot-germ-pp.R")'

# 09_plot-germ-qgparam.R ----
ms/figures/h2-germ.pdf: r/objects/fit.rds r/functions.R r/header.R r/09_plot-germ-qgparam.R
	Rscript -e 'source("r/09_plot-germ-qgparam.R")'
ms/export/vc_table_germ.rds: r/objects/fit.rds r/functions.R r/header.R r/09_plot-germ-qgparam.R
	Rscript -e 'source("r/09_plot-germ-qgparam.R")'
ms/export/diff_vpop_vg_germ.rds: r/objects/fit.rds r/functions.R r/header.R r/09_plot-germ-qgparam.R
	Rscript -e 'source("r/09_plot-germ-qgparam.R")'

# 10_plot-germ-average.R ----
r/objects/mean_germ.rds: r/objects/fit.rds r/functions.R r/header.R r/10_plot-germ-average.R
	Rscript -e 'source("r/10_plot-germ-average.R")'

# 11_plot-surv-qgparam.R ----
ms/figures/h2-surv.pdf: r/objects/fit.rds r/functions.R r/header.R r/11_plot-surv-qgparam.R
	Rscript -e 'source("r/11_plot-surv-qgparam.R")'
ms/export/vc_table_surv.rds: r/objects/fit.rds r/functions.R r/header.R r/11_plot-surv-qgparam.R
	Rscript -e 'source("r/11_plot-surv-qgparam.R")'

# 12_plot-surv-average.R ----
ms/figures/mean-traits.pdf: r/objects/fit.rds r/objects/mean_germ.rds r/functions.R r/header.R r/12_plot-surv-average.R
	Rscript -e 'source("r/12_plot-surv-average.R")'

# 13_test-selection.R ----
ms/figures/selection.pdf: processed-data/germination.rds r/objects/fit.rds r/functions.R r/header.R r/13_test-selection.R
	Rscript -e 'source("r/13_test-selection.R")'
ms/export/df_ind_lm.rds: processed-data/germination.rds r/objects/fit.rds r/functions.R r/header.R r/13_test-selection.R
	Rscript -e 'source("r/13_test-selection.R")'
ms/export/df_pop_lm.rds: processed-data/germination.rds r/objects/fit.rds r/functions.R r/header.R r/13_test-selection.R
	Rscript -e 'source("r/13_test-selection.R")'

# 14_plot-climate.R ----
ms/figures/climate.pdf: raw-data/climate_data.csv r/functions.R r/header.R r/14_plot-climate.R
	Rscript -e 'source("r/14_plot-climate.R")'

# 15_plot-range-map.R ----
ms/figures/range-map.pdf: r/functions.R r/header.R r/15_plot-range-map.R
	Rscript -e 'source("r/15_plot-range-map.R")'

# 16_archive-data.R ----
ms/muir-etal-2022.txt: r/header.R processed-data/germination.rds r/16_archive-data.R
	Rscript -e 'source("r/16_archive-data.R")'

# paper ----
ms/ms.pdf: ms/ms.Rmd ms/mimulus-germination.bib ms/export/data.rds ms/export/df_ind_lm.rds ms/export/df_pop_lm.rds ms/export/diff_vpop_vg_germ.rds ms/export/nRemove.rds ms/export/vc_table_germ.rds ms/export/vc_table_surv.rds ms/export/diff_vpop_vg_germ.rds ms/figures/climate.pdf ms/figures/h2-germ.pdf ms/figures/h2-surv.pdf ms/figures/mean-traits.pdf ms/figures/pp_check_germ.pdf ms/figures/range-map.pdf ms/figures/selection.pdf
	Rscript -e 'rmarkdown::render("ms/ms.Rmd", output_format = "bookdown::pdf_document2", output_file = "ms.pdf")'
ms/tables-and-captions.pdf: ms/ms.Rmd ms/tables-and-captions.Rmd ms/mimulus-germination.bib ms/export/data.rds ms/export/df_ind_lm.rds ms/export/df_pop_lm.rds ms/export/diff_vpop_vg_germ.rds ms/export/nRemove.rds ms/export/vc_table_germ.rds ms/export/vc_table_surv.rds ms/export/diff_vpop_vg_germ.rds ms/figures/climate.pdf ms/figures/h2-germ.pdf ms/figures/h2-surv.pdf ms/figures/mean-traits.pdf ms/figures/pp_check_germ.pdf ms/figures/range-map.pdf ms/figures/selection.pdf
	Rscript -e 'rmarkdown::render("ms/tables-and-captions.Rmd", output_format = "bookdown::pdf_document2", output_file = "tables-and-captions.pdf")'
ms/si.pdf: ms/si.Rmd ms/mimulus-germination.bib ms/export/data.rds ms/export/df_ind_lm.rds ms/export/df_pop_lm.rds ms/export/diff_vpop_vg_germ.rds ms/export/nRemove.rds ms/export/vc_table_germ.rds ms/export/vc_table_surv.rds ms/export/diff_vpop_vg_germ.rds ms/figures/climate.pdf ms/figures/h2-germ.pdf ms/figures/h2-surv.pdf ms/figures/mean-traits.pdf ms/figures/pp_check_germ.pdf ms/figures/range-map.pdf ms/figures/selection.pdf
	Rscript -e 'rmarkdown::render("ms/si.Rmd", output_format = "bookdown::pdf_document2", output_file = "si.pdf")'

clean: 
	\rm -f *~ *.Rout */*~ */*.Rout .RData Rplots.pdf
	
cleanall: 
	\rm -f *.aux *.bbl *.blg *.log *.pdf *~ *.Rout */*~ */*.Rout ms/export/*.rds ms/figures/*.pdf ms/ms.pdf ms/ms.tex ms/si.pdf ms/si.tex r/objects/*.rds processed-data/*.rds */*.aux */*.log 

