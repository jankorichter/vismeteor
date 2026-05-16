PKG     := $(shell awk '/^Package:/ {print $$2}' DESCRIPTION)
VERSION := $(shell awk '/^Version:/ {print $$2}' DESCRIPTION)
TARBALL := $(PKG)_$(VERSION).tar.gz

.DEFAULT_GOAL := build
.PHONY: styler lintr docs build check clean

styler:
	Rscript -e 'styler::style_pkg(indent_by = 4L, scope = "tokens")'
	Rscript -e 'styler::style_dir("vignettes", filetype = c("R","Rmd"), indent_by = 4L)'

lintr:
	Rscript -e 'l <- lintr::lint_package(); if (length(l) > 0) { print(l); quit(status = 1) }'

docs:
	Rscript -e 'roxygen2::roxygenise()'

build: docs lintr
	R CMD build .

check: build
	R CMD check --as-cran $(TARBALL)

clean:
	rm -rf $(PKG)_*.tar.gz $(PKG).Rcheck
