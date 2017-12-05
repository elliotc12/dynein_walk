CPPFLAGS = -std=c++98 -g -Werror -O2 -Wall
LIBRARIES = -lm

PAPERS = papers/fft_speedup/convolution_theorem.pdf papers/paper/paper.pdf papers/notes/notes.pdf papers/elliott-thesis/latex/capek.pdf
DRAW = scripts/dynein/draw/motor_domain.py scripts/dynein/draw/tail.py
HEADERS = $(wildcard src/*.h) src/version-info.h

THESIS-PLOTS = plots/trajectory-plot_thesis.pdf plots/stepping_time_histogram_thesis.pdf plots/stepping_length_histogram_thesis.pdf

PAPER-PLOTS = plots/paper-trajectory-plot.pdf plots/stepping_time_histogram_paper.pdf plots/stepping_length_histogram_paper.pdf plots/total_step_time_vs_step_length_multiple_seeds.pdf

all: generate_stepping_data public $(DRAW)

.PHONY: clean public

######### SRC stuff ##########
src/version-info.h: SHELL:=/bin/bash
src/version-info.h: .git/refs/heads/master $(wildcard src/*.cpp) Makefile
	UNAMESTR=$$(uname); if [[ "$$UNAMESTR" == 'Linux' ]]; then \
	echo -n static const char '*version' = '"' > src/version-info.h~; \
	git describe --dirty --tags | tr -d '\n' >> src/version-info.h~; \
	echo '";' >> src/version-info.h~; \
	if cmp src/version-info.h src/version-info.h~; then \
		rm src/version-info.h~; \
	else \
		echo new version needed; \
		diff -u src/version-info.h src/version-info.h~; \
		mv src/version-info.h~ src/version-info.h; \
	fi; \
	elif [[ "$$UNAMESTR" == 'Darwin' ]]; then \
	    echo 'static const char *version = "mac-breaks-version-stuff";' > src/version-info.h; \
	fi;

build/%.o: src/%.cpp $(HEADERS)
	$(CXX) -c $^ $(CPPFLAGS)
	mkdir -p build
	mv $*.o build/

generate_stepping_data: build/generate_stepping_data.o build/dynein_simulate.o \
			build/dynein_struct_onebound.o build/dynein_struct_bothbound.o \
			build/utilities.o src/version-info.h
	$(CXX) -o generate_stepping_data $^

######### draw module stuff ##########
scripts/dynein/draw/motor_domain.py: scripts/dynein/draw/create_MD_array.py scripts/dynein/draw/outer_coords.txt
	cd scripts/dynein/draw && python2 create_MD_array.py

scripts/dynein/draw/tail.py: scripts/dynein/draw/tailDomain.py
	cd scripts/dynein/draw && python2 tailDomain.py

######### data ##########
data/thesis_stepping_data.txt data/thesis_movie_data.txt: scripts/dynein/run.py scripts/generate-thesis-data.py
	python3 scripts/generate-thesis-data.py

# Taken out of make, added data file to repository:
# data/paper_trajectory_stepping_data.txt data/paper_trajectory_movie_data.txt: generate_stepping_data scripts/dynein/run.py scripts/generate-paper-trajectory-data.py
# 	python3 scripts/generate-paper-trajectory-data.py

# Taken out of make, added data file to repository:
# data/paper_histogram_stepping_data.txt: generate_stepping_data scripts/dynein/run.py scripts/histogram-helper.py
# 	python3 scripts/histogram-helper.py

######### plots ##########
plots/trajectory-plot_thesis.pdf: data/thesis_movie_data.txt scripts/trajectory-plt.py $(DRAW)
	python3 scripts/trajectory-plt.py data/thesis_movie_data.txt
	mv plots/trajectory-plot.pdf plots/trajectory-plot_thesis.pdf

plots/stepping_time_histogram_thesis.pdf plots/stepping_length_histogram_thesis.pdf: scripts/make_stepping_plots.py data/thesis_stepping_data.txt
	python3 scripts/make_stepping_plots.py data/thesis_stepping_data.txt
	mv -u plots/stepping_length_histogram.pdf plots/stepping_length_histogram_thesis.pdf
	mv -u plots/stepping_time_histogram.pdf plots/stepping_time_histogram_thesis.pdf

HISTOGRAM_DATA = $(wildcard data/paper_histogram_stepping_data*.txt)
plots/stepping_time_histogram_paper.pdf plots/stepping_length_histogram_paper.pdf: scripts/paper-histogram-plt.py $(HISTOGRAM_DATA)
	python3 scripts/paper-histogram-plt.py
	mv -u plots/stepping_length_histogram.pdf plots/stepping_length_histogram_paper.pdf
	mv -u plots/stepping_time_histogram.pdf plots/stepping_time_histogram_paper.pdf

plots/stepping_time_histogram_%.pdf plots/stepping_length_histogram_%.pdf: scripts/make_stepping_plots.py $(HISTOGRAM_DATA)
	python3 scripts/make_stepping_plots.py $*
	mv -u plots/stepping_length_histogram.pdf plots/stepping_length_histogram_$*.pdf
	mv -u plots/stepping_time_histogram.pdf plots/stepping_time_histogram_$*.pdf

plots/total_step_time_vs_step_length_multiple_seeds.pdf: scripts/color_hist.py $(HISTOGRAM_DATA)
	python scripts/color_hist.py -v -a

plots/paper-trajectory-plot.pdf: data/paper_trajectory_movie_data.txt scripts/paper-trajectory-plt.py $(DRAW)
	python3 scripts/paper-trajectory-plt.py data/paper_trajectory

plots/paper-movie.mp4: data/paper_trajectory_movie_data.txt scripts/movie.py $(DRAW)
	mkdir -p movies
	./scripts/movie.py data/paper_trajectory_movie_data.txt speed=1 tail forces
	mv plots/movie.mp4 plots/paper-movie.mp4

######### papers ##########
SVG_FIGURES = $(wildcard papers/*/figures/*.svg)
FIGURES = $(patsubst %.svg,%.pdf,$(SVG_FIGURES))
ALL_FIGURES = $(FIGURES) $(THESIS-PLOTS) $(PAPER-PLOTS)

papers/elliott-thesis/figures/%.pdf: papers/elliott-thesis/figures/%.svg
	inkscape -D --export-pdf $(shell pwd)/$@ $(shell pwd)/$<

papers/paper/figures/%.pdf: papers/paper/figures/%.svg
	inkscape -D --export-pdf $(shell pwd)/$@ $(shell pwd)/$<

papers/fft_speedup/convolution_theorem.pdf: papers/fft_speedup/convolution_theorem.tex
	cd papers/fft_speedup && pdflatex convolution_theorem &&  pdflatex convolution_theorem

papers/elliott-thesis/latex/capek.pdf: papers/elliott-thesis/latex/thesis.tex $(ALL_FIGURES)
	cd papers/elliott-thesis/latex && xelatex thesis.tex && bibtex thesis && xelatex thesis.tex && xelatex thesis.tex
	mv papers/elliott-thesis/latex/thesis.pdf papers/elliott-thesis/latex/capek.pdf

papers/paper/paper.pdf: papers/paper/paper.tex $(ALL_FIGURES)
	(cd papers/paper && xelatex paper.tex && bibtex paper && xelatex paper.tex && xelatex paper.tex) || (rm -f $@ && false)

papers/notes/notes.pdf: papers/notes/notes.tex
	cd papers/notes && xelatex notes.tex

public: $(PAPERS)
	cp -v $(PAPERS) public/

clean:
	rm -f build/*.o
	rm -f src/version-info.h
	rm -f generate_stepping_data
	rm -f scripts/dynein/draw/motor_domain.py scripts/dynein/draw/tail.py
	rm -f scripts/*.pyc scripts/*/*.pyc
	rm -rf plots
	rm -f $(ALL_FIGURES)
	rm -f public/*.pdf public/*.tex
	rm -f data/thesis_*
