# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER 		  =
SPHINXPROJ    = opcsim
SOURCEDIR     = .
BUILDDIR      = _build

# Internal Variables
PAPEROPT_a4 	= -D latex_paper_size=a4
PAPEROPT_letter	= -D latex_paper_size=letter
ALLSPHINXOPTS 	= -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help clean html Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	@echo "Cleaning..."
	-rm -rf $(BUILDDIR)/*
	-rm -rf examples/*
	-rm -rf example_thumbs/*
	-rm -rf tutorial/*_files/
	-rm -rf tutorial/*.rst
	-rm -rf generated/*

tutorials:
	make -C tutorial

introduction:
	tools/nb_to_doc.py introduction

notebooks: tutorials introduction

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."


upload:
	@echo "Uploading to github@gh-pages"
