CPPFLAGS = -std=c++98 -g -Werror -O2 -Wall
LIBRARIES = -lm

PAPERS = papers/fft_speedup/convolution_theorem.pdf
DRAW = scripts/draw/motor_domain.py scripts/draw/tail.py
HEADERS = $(wildcard src/*.h) src/version-info.h

all: generate_stepping_data $(PAPERS) $(DRAW)

.PHONY: clean elliott-thesis

.PRECIOUS: data/stepping_data_%.txt data/stepping_config_%.txt data/stepping_movie_data_%.txt data/bothbound_data_%.bin data/onebound_data_%.bin data/ob_config_%.txt data/bb_config_%.txt data/exploration_movie_data.txt data/exploration_stepping_data.txt # prevent nonexistant data files from being deleted after creation+use


######### SRC stuff ##########
src/version-info.h: SHELL:=/bin/bash
src/version-info.h: ../.git/refs/heads/master $(wildcard src/*.cpp) Makefile
	UNAMESTR=$$(uname); if [[ "$$UNAMESTR" == 'Linux' ]]; then \
	echo -n static const char '*version' = '"' > src/version-info.h; \
	git describe --dirty --tags | tr -d '\n' >> src/version-info.h; \
	echo -n '-=-' >> src/version-info.h; \
	date -Ins | tr -d '\n' >> src/version-info.h; \
	echo '";' >> src/version-info.h; \
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
scripts/draw/motor_domain.py: scripts/draw/create_MD_array.py scripts/draw/outer_coords.txt
	cd scripts/draw && python2 create_MD_array.py

scripts/draw/tail.py: scripts/draw/tailDomain.py
	cd scripts/draw && python2 tailDomain.py

######### data ##########
data/thesis_stepping_data.txt data/thesis_movie_data.txt: generate_stepping_data scripts/run_scripts/simrunner.py scripts/run_scripts/generate-thesis-data.py
	mkdir -p data
	mkdir -p scripts/run_scripts/runlogs
	python3 scripts/run_scripts/generate-thesis-data.py

data/exploration_stepping_data.txt data/exploration_movie_data.txt: generate_stepping_data scripts/run_scripts/simrunner.py scripts/run_scripts/generate-exploration-data.py
	mkdir -p data
	mkdir -p scripts/run_scripts/runlogs
	python3 scripts/run_scripts/generate-exploration-data.py

######### plots ##########
plots/trajectory-plot_thesis.pdf: data/thesis_movie_data.txt scripts/trajectory-plt.py $(wildcard scripts/draw/*.py)
	python3 scripts/trajectory-plt.py data/thesis_movie_data.txt
	mv plots/trajectory-plot.pdf plots/trajectory-plot_thesis.pdf

plots/stepping_time_histogram_thesis.pdf plots/stepping_length_histogram_thesis.pdf: scripts/make_stepping_plots.py data/thesis_stepping_data.txt
	python3 scripts/make_stepping_plots.py data/thesis_stepping_data.txt
	mv plots/stepping_length_histogram.pdf plots/stepping_length_histogram_thesis.pdf
	mv plots/stepping_time_histogram.pdf plots/stepping_time_histogram_thesis.pdf

plots/exploration-plot.pdf plots/stepping_time_histogram_exploration.pdf plots/stepping_length_histogram_exploration.pdf: data/exploration_movie_data.txt scripts/exploration-plt.py $(DRAW)
	python3 scripts/exploration-plt.py data/exploration_movie_data.txt data/exploration_stepping_data.txt
	mv plots/stepping_length_histogram.pdf plots/stepping_length_histogram_exploration.pdf
	mv plots/stepping_time_histogram.pdf plots/stepping_time_histogram_exploration.pdf

plots/exploration-movie.mp4: data/exploration_movie_data.txt scripts/movie.py $(DRAW)
	mkdir -p movies
	./scripts/movie.py data/exploration_movie_data.txt speed=1 tail
	mv plots/movie.mp4 plots/exploration-movie.mp4

######### papers ##########
SVG_FIGURES = $(wildcard papers/*/figures/*.svg)

ALL_FIGURES = $(patsubst %.svg,%.pdf,$(SVG_FIGURES))

papers/elliott-thesis/figures/%.pdf: papers/elliott-thesis/figures/%.svg
	inkscape -D --export-pdf $(shell pwd)/$@ $(shell pwd)/$<

papers/fft_speedup/convolution_theorem.pdf: papers/fft_speedup/convolution_theorem.tex
	cd papers/fft_speedup && pdflatex convolution_theorem &&  pdflatex convolution_theorem

elliott-thesis: $(DRAW) $(ALL_FIGURES)
	$(MAKE) -C papers/elliott-thesis

clean:
	rm -f build/*.o
	rm -f src/version-info.h
	rm -f generate_stepping_data
	rm -f scripts/draw/motor_domain.py scripts/draw/tail.py
	rm -f scripts/*.pyc scripts/*/*.pyc
	rm -rf plots
	$(MAKE) -C papers/elliott-thesis clean
	rm -f $(patsubst %.svg,%.pdf,$(SVG_FIGURES))
