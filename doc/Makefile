
doc = sxamg.pdf

SRC_TEX   = $(wildcard *.tex)
TEXFILES = $(SRC_TEX)

LATEX = latex
DVIPDF = dvipdf

default: $(doc)

sxamg.pdf: sxamg.dvi
	$(DVIPDF) sxamg.dvi

sxamg.dvi: $(TEXFILES)
	$(LATEX) sxamg.tex
	$(LATEX) sxamg.tex

.PHONY: clean

clean:
	@rm -fv *.aux sxamg.dvi *.log *.out *.toc sxamg.idx  sxamg.ind sxamg.ilg
