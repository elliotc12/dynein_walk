FIGURES = $(wildcard *.svg)

TO_DELETE=$(patsubst %.svg, %.pdf, $(wildcard *.svg))

all: $(patsubst %.svg,%.pdf,$(FIGURES))

%.pdf: %.svg
	inkscape -D --export-pdf $@ $<
# doesn't work, fix sometime

clean:
	rm -f $(TO_DELETE)
