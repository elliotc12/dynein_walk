FIGURES = $(wildcard *.svg)

all: $(patsubst %.svg,%.pdf,$(FIGURES))

%.pdf: %.svg
	inkscape -D --export-pdf $@ $<
# doesn't work, fix sometime
