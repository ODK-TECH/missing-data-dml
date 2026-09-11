# Reproducible pipeline. Run `make all` from the repository root.

RSCRIPT = Rscript

.PHONY: all deps smoke sim report clean

all: sim report

deps:
	$(RSCRIPT) scripts/install-deps.R

smoke:
	$(RSCRIPT) tests/test-smoke.R

sim:
	$(RSCRIPT) scripts/run-all.R

quick:
	QUICK=1 $(RSCRIPT) scripts/run-all.R

report:
	$(RSCRIPT) -e 'rmarkdown::render("report/report.Rmd", output_format = "html_document")'

clean:
	rm -f results/*.rds results/*.csv results/session-info.txt figures/*.png report/report.html
