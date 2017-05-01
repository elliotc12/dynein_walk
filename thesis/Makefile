default: thesis.pdf

FIGURES=$(patsubst ../figures/%.svg, ../figures/%.pdf, $(wildcard ../figures/*.svg))

PLOTS= ../plots/x-trajectory_thesis.pdf ../plots/y-trajectory_thesis.pdf ../plots/stepping_length_histogram_thesis.pdf ../plots/stepping_time_histogram_thesis.pdf

thesis.pdf: latex/thesis.tex $(FIGURES) $(PLOTS)
	#cd latex && xelatex thesis.tex && bibtex thesis && xelatex thesis.tex && xelatex thesis.tex && mv thesis.pdf ..
	cd latex && xelatex -interaction=batchmode thesis.tex && tail -n 25 thesis.log && mv thesis.pdf ..

# thesis.pdf: latex/thesis.tex $(FIGURES)
# 	cd latex && xelatex thesis.tex && mv thesis.pdf ..

thesis2.pdf: latex/thesis2.tex $(FIGURES)
	cd latex && xelatex thesis2.tex && mv thesis2.pdf ..

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
