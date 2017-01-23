default: thesis.pdf

FIGURES=$(patsubst ../figures/%.svg, ../figures/%.pdf, $(wildcard ../figures/*.svg))

thesis.pdf: latex/thesis.tex $(FIGURES)
	cd latex && xelatex thesis.tex && bibtex thesis && xelatex thesis.tex && xelatex thesis.tex && mv thesis.pdf ..

continuous:
	touch latex/thesis.pdf
	ln -sf latex/thesis.pdf thesis.pdf
	cd latex && latexmk -pvc --xelatex thesis.tex

../figures/%.pdf: ../figures/%.svg
	cd ../figures && make

clean:
	rm -f latex/*.fdb_latexmk
	rm -f latex/*.log
	rm -f latex/*.pdf
	rm -f latex/*.aux
	rm -f latex/*.fls
	rm -f latex/*.out
	rm -f latex/*.toc
	rm -f latex/*.blg
	rm -f latex/*.bbl
	rm -f *.pdf
