CPPFLAGS = -std=c++98 -g -Werror -O2 -Wall
LIBRARIES = -lm

PAPERS = papers/fft_speedup/convolution_theorem.pdf papers/paper/paper.pdf papers/notes/notes.pdf
DRAW = scripts/dynein/draw/motor_domain.py scripts/dynein/draw/tail.py
HEADERS = $(wildcard src/*.h)

PARAMSEARCH_DATAFILES = $(wildcard data/parameterSearch/*.txt)
PARAMSEARCH_PDFS = $(patsubst data/parameterSearch/%.txt,plots/parameterSearch/%.pdf,$(PARAMSEARCH_DATAFILES))

RANDOMSEARCH_DATAFILES = $(wildcard data/stepping_data_randomsearch__*.txt)
RANDOMSEARCH_PDFS = $(patsubst data/stepping_data_randomsearch__%.txt, plots/randomsearch/%.pdf,$(RANDOMSEARCH_DATAFILES))

UNBINDING_PROBABILITY_PDFS = $(patsubst data/unbinding_probability/%.tex,plots/unbinding_probability/%.pdf, $(wildcard data/unbinding_probability/*.tex))

all: generate_stepping_data public $(DRAW) $(PAPERS) $(RANDOMSEARCH_PDFS) # $(PARAMSEARCH_PDFS)

.PHONY: clean public

######### SRC stuff ##########
build/%.o: src/%.cpp $(HEADERS)
	$(CXX) -c $^ $(CPPFLAGS)
	mkdir -p build
	mv $*.o build/

generate_stepping_data: build/generate_stepping_data.o build/dynein_simulate.o \
			build/dynein_struct_onebound.o build/dynein_struct_bothbound.o \
			build/utilities.o
	$(CXX) -o generate_stepping_data $^

simulate_unbinding_rates: build/simulate_unbinding_rates.o build/dynein_simulate.o \
			build/dynein_struct_onebound.o build/dynein_struct_bothbound.o \
			build/utilities.o
	$(CXX) -o simulate_unbinding_rates $^

######### draw module stuff ##########
scripts/dynein/draw/motor_domain.py: scripts/dynein/draw/create_MD_array.py scripts/dynein/draw/outer_coords.txt
	cd scripts/dynein/draw && python create_MD_array.py

scripts/dynein/draw/tail.py: scripts/dynein/draw/tailDomain.py
	cd scripts/dynein/draw && python tailDomain.py

######### data ##########
data/thesis_stepping_data.txt data/thesis_movie_data.txt: scripts/dynein/run.py scripts/generate-thesis-data.py
	python3 scripts/generate-thesis-data.py

data/paper_trajectory_stepping_data.txt data/paper_trajectory_movie_data.txt: generate_stepping_data scripts/dynein/run.py scripts/generate-paper-trajectory-data.py data/paper_params.py
	python3 scripts/generate-paper-trajectory-data.py

data/paper_unbinding_probability__L-10,s-1.txt: simulate_unbinding_rates scripts/generate-paper-unbinding-rate-data.py data/paper_params.py
	rm -f data/paper_unbinding_probability__*
	python3 scripts/generate-paper-unbinding-rate-data.py

# Taken out of make, added data file to repository:
# data/paper_histogram_stepping_data.txt: generate_stepping_data scripts/dynein/run.py scripts/histogram-helper.py
# 	python3 scripts/histogram-helper.py

######### fun plots ##########
plots/stepping_time_histogram_paper.pdf plots/stepping_length_histogram_paper.pdf plots/stepping_analysis_paper.pdf: scripts/paper-histogram-plt.py $(STATIC_DATA)
	python3 scripts/paper-histogram-plt.py
	mv plots/stepping_length_histogram.pdf plots/stepping_length_histogram_paper.pdf
	mv plots/stepping_time_histogram.pdf plots/stepping_time_histogram_paper.pdf
	mv plots/stepping_analysis.pdf plots/stepping_analysis_paper.pdf

plots/stepping_time_histogram_%.pdf: scripts/make_stepping_plots.py
	python3 scripts/make_stepping_plots.py $*
	mv plots/stepping_time_histogram.pdf plots/stepping_time_histogram_$*.pdf

######### paper plots ##########
PAPER_DATA = $(wildcard data/paper_main_stepping_data*.txt)

PAPER-PLOTS = plots/paper_trajectory_plot.pdf plots/paper_displacement_vs_step_length.pdf plots/paper_onebound_vs_steplength.pdf plots/paper_initial_vs_final_displacement.pdf plots/paper_foot_order_histogram.pdf plots/paper-trajectory-movie.mp4 plots/paper_model_behavior.pdf

plots/paper_model_behavior.pdf plots/paper_foot_order_histogram.pdf: scripts/make_all_stepping_plots.py $(PAPER_DATA) data/paper_unbinding_probability__L-10,s-1.txt
	python3 scripts/make_all_stepping_plots.py -d data -b paper_main
	mv plots/stepping_analysis.pdf plots/paper_foot_order_histogram.pdf
	mv plots/model_behavior.pdf plots/paper_model_behavior.pdf

plots/paper_initial_vs_final_displacement.pdf plots/paper_onebound_vs_steplength.pdf plots/paper_displacement_vs_step_length.pdf: scripts/color_hist2.py $(PAPER_DATA)
	python3 scripts/color_hist2.py -d paper_main
	mv plots/initial-vs-final-displacement.pdf plots/paper_initial_vs_final_displacement.pdf
	mv plots/onebound-time-vs-step-length.pdf plots/paper_onebound_vs_steplength.pdf
	mv plots/initial-displacement-vs-step-length.pdf plots/paper_displacement_vs_step_length.pdf

plots/paper_trajectory_plot.pdf: data/paper_trajectory_movie_data.txt scripts/paper-trajectory-plt.py $(DRAW)
	python3 scripts/paper-trajectory-plt.py data/paper_trajectory

plots/burgess-model-figure.pdf plots/grotjahn-model-figure.pdf: scripts/generate-paper-model-figures.py papers/paper/figures/model-raw-images/burgess-fig-4-cropped.png papers/paper/figures/model-raw-images/grotjahn-model-figure.png data/paper_params.tex
	python3 scripts/generate-paper-model-figures.py

data/paper_params.tex: data/paper_params.py scripts/parameters-to-latex.py
	python3 scripts/parameters-to-latex.py

plots/paper-trajectory-movie.mp4: data/paper_trajectory_movie_data.txt scripts/movie.py
	python2.7 scripts/movie.py data/paper_trajectory_movie_data.txt speed=1000
	mv plots/movie.mp4 plots/paper-trajectory-movie.mp4

######### thesis plots ##########
THESIS-PLOTS = plots/trajectory-plot_thesis.pdf plots/stepping_time_histogram_thesis.pdf plots/stepping_length_histogram_thesis.pdf

plots/trajectory-plot_thesis.pdf: data/thesis_movie_data.txt scripts/trajectory-plt.py $(DRAW)
	python3 scripts/trajectory-plt.py data/thesis_movie_data.txt
	mv plots/trajectory-plot.pdf plots/trajectory-plot_thesis.pdf

plots/stepping_time_histogram_thesis.pdf plots/stepping_length_histogram_thesis.pdf: scripts/make_thesis_stepping_plots.py data/thesis_stepping_data.txt
	python3 scripts/make_thesis_stepping_plots.py data/thesis_stepping_data.txt
	mv plots/stepping_length_histogram.pdf plots/stepping_length_histogram_thesis.pdf
	mv plots/stepping_time_histogram.pdf plots/stepping_time_histogram_thesis.pdf

######### parameterSearch PDFs ##########
plots/parameterSearch/%.pdf: data/parameterSearch/%.txt data/parameterSearch/%.tex scripts/make_all_stepping_plots.py scripts/color_hist2.py plots/parameterSearch/display_template.tex
	mkdir -p plots/parameterSearch/searchplots
	python3 scripts/make_all_stepping_plots.py -d data/parameterSearch -b $*
	python3 scripts/color_hist2.py -r data/parameterSearch/ -d $*
	cp data/parameterSearch/$*.tex plots/parameterSearch/search_parameters.tex
	cd plots/parameterSearch && xelatex display_template.tex
	mv plots/parameterSearch/display_template.pdf plots/parameterSearch/$*.pdf

plots/randomsearch/%.pdf: data/stepping_data_randomsearch__%.txt scripts/make_all_stepping_plots.py scripts/color_hist2.py plots/parameterSearch/display_template.tex
	mkdir -p plots/randomsearch/
	python3 scripts/make_all_stepping_plots.py -q -d data -b $*
	mv plots/model_behavior.pdf plots/randomsearch/$*.pdf

######### unbinding probability PDFs ##########
plots/unbinding_probability/%.pdf: $(wildcard data/unbinding_probability/%*) data/unbinding_probability/%.tex scripts/make_all_unbinding_probability_plots.py plots/unbinding_probability/display_template.tex
	mkdir -p plots/unbinding_probability/plots_for_latex
	python3 scripts/make_all_unbinding_probability_plots.py -d data/unbinding_probability -b $*
	cp data/unbinding_probability/$*.tex plots/unbinding_probability/parameters.tex
	cd plots/unbinding_probability && xelatex display_template.tex
	mv plots/unbinding_probability/display_template.pdf plots/unbinding_probability/$*.pdf

######### papers ##########
PAPER_SVG_FIGURES = $(wildcard papers/*/figures/*.svg)
PAPER-FIGURES = $(patsubst %.svg,%.pdf,$(PAPER_SVG_FIGURES)) plots/burgess-model-figure.pdf plots/grotjahn-model-figure.pdf

papers/elliott-thesis/figures/%.pdf: papers/elliott-thesis/figures/%.svg
	inkscape -D --export-pdf $(shell pwd)/$@ $(shell pwd)/$<

papers/paper/figures/%.pdf: papers/paper/figures/%.svg
	inkscape -D --export-pdf $(shell pwd)/$@ $(shell pwd)/$<

papers/fft_speedup/convolution_theorem.pdf: papers/fft_speedup/convolution_theorem.tex
	cd papers/fft_speedup && pdflatex convolution_theorem &&  pdflatex convolution_theorem

papers/elliott-thesis/latex/capek.pdf: papers/elliott-thesis/latex/thesis.tex $(THESIS-PLOTS)
	cd papers/elliott-thesis/latex && xelatex thesis.tex && bibtex thesis && xelatex thesis.tex && xelatex thesis.tex
	mv papers/elliott-thesis/latex/thesis.pdf papers/elliott-thesis/latex/capek.pdf

papers/paper/paper.pdf: papers/paper/paper.tex papers/paper/paper-supplemental.tex $(PAPER-FIGURES) $(PAPER-PLOTS) data/paper_params.tex
	(cd papers/paper && xelatex paper.tex && bibtex paper && xelatex paper.tex && xelatex paper.tex) || (rm -f $@ && false)
	(cd papers/paper && xelatex paper-supplemental.tex && xelatex paper-supplemental.tex) || (rm -f $@ && false)

papers/notes/notes.pdf: papers/notes/notes.tex
	cd papers/notes && xelatex notes.tex

public: $(PAPERS)
	cp -v $(PAPERS) public/

clean:
	rm -f build/*.o
	rm -f generate_stepping_data
	rm -f scripts/dynein/draw/motor_domain.py scripts/dynein/draw/tail.py
	rm -f scripts/*.pyc scripts/*/*.pyc
	rm -rf plots
	rm -f $(PAPER-FIGURES) $(THESIS-PLOTS) $(PAPER-PLOTS)
	rm -f public/*.pdf public/*.tex
	rm -f data/thesis_*
