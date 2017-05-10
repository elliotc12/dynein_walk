CPPFLAGS = -std=c++98 -g -Werror -O2 -Wall
LIBRARIES = -lm

FIGURES=$(patsubst %.svg,%.pdf,$(wildcard figures/*.svg))
STEPPING_LENGTH_HISTOGRAMS=$(patsubst data/stepping_data_%.txt, plots/stepping_length_histogram_%.pdf, $(wildcard data/stepping_data_*.txt))
STEPPING_TIME_HISTOGRAMS=$(patsubst data/stepping_data_%.txt, plots/stepping_time_histogram_%.pdf, $(wildcard data/stepping_data_*.txt))

STEPPING_MOVIES=$(patsubst data/stepping_movie_data_%.txt, movies/%.mp4, $(wildcard data/stepping_movie_data_*.txt))

.PHONY: clean histograms thesis

.PRECIOUS: data/stepping_data_%.txt data/stepping_config_%.txt data/stepping_movie_data_%.txt data/bothbound_data_%.bin data/onebound_data_%.bin data/ob_config_%.txt data/bb_config_%.txt # prevent nonexistant data files from being deleted after creation+use

all: test_onebound.results test_bothbound.results create_ob_plots create_ob_movie thesis_stuff.pdf generate_stepping_data thesis

histogram-stuff: test_onebound.results test_bothbound.results generate_stepping_data # called by make_histograms.py

version-info.h: SHELL:=/bin/bash
version-info.h: .git/refs/heads/master $(wildcard simulation/*.cpp) Makefile
	UNAMESTR=$$(uname); if [[ "$$UNAMESTR" == 'Linux' ]]; then \
	echo -n static const char '*version' = '"' > version-info.h; \
	cat version-info.h; \
	git describe --dirty --tags | tr -d '\n' >> version-info.h; \
	cat version-info.h; \
	echo -n '-=-' >> version-info.h; \
	cat version-info.h; \
	date -Ins | tr -d '\n' >> version-info.h; \
	cat version-info.h; \
	echo '";' >> version-info.h; \
	cat version-info.h; \
	elif [[ "$$UNAMESTR" == 'Darwin' ]]; then \
	    echo 'static const char *version = "mac-breaks-version-stuff";' > version-info.h; \
	fi;

derivation.pdf: latex/derivation.tex $(FIGURES)
	cd latex && pdflatex derivation.tex && mv derivation.pdf ..

derivation_confirmation.pdf: latex/derivation_confirmation.tex $(FIGURES)
	cd latex && pdflatex derivation_confirmation.tex && mv derivation_confirmation.pdf ..

test_bothbound.o: test_bothbound.cpp dynein_struct.h default_parameters.h
	$(CXX) -c test_bothbound.cpp $(CPPFLAGS)

dynein_struct_onebound.o: dynein_struct_onebound.cpp dynein_struct.h default_parameters.h
	$(CXX) -c dynein_struct_onebound.cpp $(CPPFLAGS)

dynein_struct_bothbound.o: dynein_struct_bothbound.cpp dynein_struct.h default_parameters.h
	$(CXX) -c dynein_struct_bothbound.cpp $(CPPFLAGS)

dynein_simulate.o: dynein_simulate.cpp dynein_struct_onebound.cpp dynein_struct_bothbound.cpp default_parameters.h simulations/simulation_defaults.h dynein_struct.h
	$(CXX) -c dynein_simulate.cpp $(CPPFLAGS)

simulations.o: simulations/simulations.cpp dynein_struct.h default_parameters.h
	$(CXX) -c simulations/simulations.cpp $(CPPFLAGS) -o simulations.o

test_onebound.o: test_onebound.cpp dynein_struct.h default_parameters.h
	$(CXX) -c test_onebound.cpp $(CPPFLAGS)

figures/%.pdf: figures/Makefile figures/%.svg
	rm -f $@
	cd figures && $(MAKE) $(patsubst figures/%,%,$@)

paper.pdf: latex/paper.tex $(FIGURES)
	cd latex && pdflatex -interaction nonstopmode -halt-on-error paper.tex && mv paper.pdf ..

test_bothbound: test_bothbound.o dynein_struct_bothbound.o dynein_struct_onebound.o utilities.o dynein_simulate.o
	$(CXX) test_bothbound.o dynein_struct_bothbound.o dynein_struct_onebound.o dynein_simulate.o utilities.o -o test_bothbound

test_onebound: test_onebound.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o dynein_simulate.o
	$(CXX) test_onebound.o dynein_struct_onebound.o dynein_struct_bothbound.o dynein_simulate.o utilities.o -o test_onebound

test_bothbound.results: test_bothbound
	./test_bothbound > test_bothbound.results.failed
	mv test_bothbound.results.failed test_bothbound.results

test_onebound.results: test_onebound
	./test_onebound > test_onebound.results.failed
	mv test_onebound.results.failed test_onebound.results

utilities.o: utilities.cpp dynein_struct.h default_parameters.h simulations/simulation_defaults.h simulations/plotting_defaults.h
	$(CXX) -c utilities.cpp $(CPPFLAGS)

######################### SIMULATION STUFF ###############################
TITLE = defaultplot

simulations/simulation_results/binding_time_fraction.txt: simulations/get_binding_time_fraction.py generate_stepping_data
	mkdir -p simulations/simulation_results
	./simulations/get_binding_time_fraction.py

simulations/simulation_defaults.h: simulations/custom_simulation_parameters.h

simulations/custom_simulation_parameters.h:
	touch simulations/custom_simulation_parameters.h

create_ob_plots: simulations/create_ob_plots.cpp dynein_struct.h default_parameters.h simulations/simulation_defaults.h simulations/plotting_defaults.h dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o
	$(CXX) -c simulations/create_ob_plots.cpp $(CPPFLAGS)
	$(CXX) dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o create_ob_plots.o simulations.o -o create_ob_plots

create_bb_plots: simulations/create_bb_plots.cpp dynein_struct.h default_parameters.h simulations/simulation_defaults.h simulations/plotting_defaults.h dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o
	$(CXX) -c simulations/create_bb_plots.cpp $(CPPFLAGS)
	$(CXX) dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o create_bb_plots.o simulations.o -o create_bb_plots

create_ob_movie: simulations/create_ob_movie.cpp dynein_struct.h default_parameters.h simulations/simulation_defaults.h simulations/plotting_defaults.h dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o
	$(CXX) -c simulations/create_ob_movie.cpp $(CPPFLAGS)
	$(CXX) dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o create_ob_movie.o -o create_ob_movie

create_bb_movie: simulations/create_bb_movie.cpp dynein_struct.h default_parameters.h simulations/simulation_defaults.h simulations/plotting_defaults.h dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o
	$(CXX) -c simulations/create_bb_movie.cpp $(CPPFLAGS)
	$(CXX) dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o create_bb_movie.o -o create_bb_movie

plots/OB_Force_x_%.pdf plots/OB_Force_y_%.pdf: ob_plots.sh create_ob_plots make_plot.py data/ob_config_%.txt data/onebound_data_%.bin
	sh ob_plots.sh $*

plots/BB_Force_x_%.pdf plots/BB_Force_y_%.pdf: bb_plots.sh create_bb_plots make_plot.py data/bb_config_%.txt data/bothbound_data_%.bin
	sh bb_plots.sh $*

histograms:
	make $(STEPPING_LENGTH_HISTOGRAMS)
	make $(STEPPING_TIME_HISTOGRAMS)

generate_stepping_data: simulations/generate_stepping_data.cpp dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o test_onebound.results test_bothbound.results version-info.h
	mkdir -p runlogs
	mkdir -p data
	$(CXX) -c simulations/generate_stepping_data.cpp $(CPPFLAGS)
	$(CXX) generate_stepping_data.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o -o generate_stepping_data

plots/stepping_time_histogram_thesis.pdf plots/stepping_length_histogram_thesis.pdf: make_stepping_plots.py data/thesis_stepping_data.txt
	python3 make_stepping_plots.py data/thesis_stepping_data.txt
	mv plots/stepping_length_histogram.pdf plots/stepping_length_histogram_thesis.pdf
	mv plots/stepping_time_histogram.pdf plots/stepping_time_histogram_thesis.pdf

data/thesis_stepping_data.txt data/thesis_movie_data.txt: generate_stepping_data run_scripts/simrunner.py run_scripts/generate-thesis-data.py
	python3 run_scripts/generate-thesis-data.py

plots/trajectory-plot_thesis.pdf: data/thesis_movie_data.txt trajectory-plt.py draw_cartoon.py
	python3 trajectory-plt.py data/thesis_movie_data.txt
	mv plots/trajectory-plot.pdf plots/trajectory-plot_thesis.pdf

#data/stepping_config_%.txt data/stepping_data_%.txt data/stepping_movie_data_%.txt:
#	mkdir -p data
#	make generate_stepping_data
#	./generate_stepping_data $*

movies/ob_%.gif: create_ob_movie data/onebound_data_%.bin movie.py
	@echo "Use TITLE='yourtitle' to give plot a title"
	./create_ob_movie $*
	mkdir -p movies
	./movie.py $* speed=1

movies/bb_%.gif: create_bb_movie data/bothbound_data_%.bin movie.py
	@echo "Use TITLE='yourtitle' to give plot a title"
	./create_bb_movie $*
	mkdir -p movies
	./movie.py $* speed=1

movies/%.mp4: data/stepping_movie_data_%.txt movie.py
	@echo "Use TITLE='yourtitle' to give plot a title"
	mkdir -p movies
	./movie.py $* speed=10 tail

stepping_movies:
	make $(STEPPING_MOVIES)

data/ob_config_%.txt data/onebound_data_%.bin: #dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations/generate_onebound_data.cpp default_parameters.h dynein_struct.h simulations/simulation_defaults.h
	mkdir -p data
	$(CXX) -c simulations/generate_onebound_data.cpp $(CPPFLAGS)
	$(CXX) generate_onebound_data.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o -o generate_onebound_data
	./generate_onebound_data $*

data/bb_config_%.txt data/bothbound_data_%.bin: #dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations/generate_bothbound_data.cpp default_parameters.h dynein_struct.h simulations/simulation_defaults.h
	mkdir -p data
	$(CXX) -c simulations/generate_bothbound_data.cpp $(CPPFLAGS)
	$(CXX) generate_bothbound_data.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o -o generate_bothbound_data
	./generate_bothbound_data $*

# ob_PE_equipartition_ratio_average_vs_spring_constant: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/simulation_defaults.h simulations/ob_PE_equipartition_ratio_average_vs_spring_constant.cpp
# 	$(CXX) -c simulations/ob_PE_equipartition_ratio_average_vs_spring_constant.cpp $(CPPFLAGS)
# 	$(CXX) ob_PE_equipartition_ratio_average_vs_spring_constant.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_average_vs_spring_constant

# ob_PE_equipartition_ratio_average_vs_spring_constant_plot: ob_PE_equipartition_ratio_average_vs_spring_constant
# 	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
# 	mkdir -p plots
# 	mkdir -p data
# 	./ob_PE_equipartition_ratio_average_vs_spring_constant $(TITLE)
# 	./make_plot.py --figtitle="$(TITLE)" --logx --xlabel="Spring constant (nm^2*kg/s^2)" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/bba_pe_equipartition_ratio_vs_c_$(TITLE).txt data/bma_pe_equipartition_ratio_vs_c_$(TITLE).txt data/ta_pe_equipartition_ratio_vs_c_$(TITLE).txt data/uma_pe_equipartition_ratio_vs_c_$(TITLE).txt data/config_pe_equipartition_ratio_vs_c_$(TITLE).txt

# ob_PE_equipartition_ratio_average_vs_force_ratio: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_equipartition_ratio_average_vs_force_ratio.cpp simulations/simulation_defaults.h
# 	$(CXX) -c simulations/ob_PE_equipartition_ratio_average_vs_force_ratio.cpp $(CPPFLAGS)
# 	$(CXX) ob_PE_equipartition_ratio_average_vs_force_ratio.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_average_vs_force_ratio

# ob_PE_equipartition_ratio_average_vs_force_ratio_plot: ob_PE_equipartition_ratio_average_vs_force_ratio
# 	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
# 	mkdir -p plots
# 	mkdir -p data
# 	./ob_PE_equipartition_ratio_average_vs_force_ratio $(TITLE)
# 	./make_plot.py --scatter --figtitle="$(TITLE)" --xlabel="Brownian / conformational force variance" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/bba_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt data/bma_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt data/ta_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt data/uma_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt data/config_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt

# ob_PE_equipartition_ratio_average_vs_temperature: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_equipartition_ratio_average_vs_temperature.cpp simulations/simulation_defaults.h
# 	$(CXX) -c simulations/ob_PE_equipartition_ratio_average_vs_temperature.cpp $(CPPFLAGS)
# 	$(CXX) ob_PE_equipartition_ratio_average_vs_temperature.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_average_vs_temperature

# ob_PE_equipartition_ratio_average_vs_temperature_plot: ob_PE_equipartition_ratio_average_vs_temperature
# 	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
# 	mkdir -p plots
# 	mkdir -p data
# 	./ob_PE_equipartition_ratio_average_vs_temperature $(TITLE)
# 	./make_plot.py --figtitle="$(TITLE)" --xlabel="Temp (K)" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/bba_pe_equipartition_ratio_vs_T_$(TITLE).txt data/bma_pe_equipartition_ratio_vs_T_$(TITLE).txt data/ta_pe_equipartition_ratio_vs_T_$(TITLE).txt data/uma_pe_equipartition_ratio_vs_T_$(TITLE).txt data/config_pe_equipartition_ratio_vs_T_$(TITLE).txt

# ob_PE_equipartition_ratio_average_vs_dt: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_equipartition_ratio_average_vs_dt.cpp
# 	$(CXX) -c simulations/ob_PE_equipartition_ratio_average_vs_dt.cpp $(CPPFLAGS)
# 	$(CXX) ob_PE_equipartition_ratio_average_vs_dt.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_average_vs_dt

# ob_PE_equipartition_ratio_average_vs_dt_plot: ob_PE_equipartition_ratio_average_vs_dt
# 	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
# 	mkdir -p plots
# 	mkdir -p data
# 	./ob_PE_equipartition_ratio_average_vs_dt $(TITLE)
# 	./make_plot.py --figtitle="$(TITLE)" --xlabel="log(dt)" --ylabel="PE / 0.5*kb*T" --logx data/bba_pe_equipartition_ratio_vs_dt_$(TITLE).txt data/bma_pe_equipartition_ratio_vs_dt_$(TITLE).txt data/ta_pe_equipartition_ratio_vs_dt_$(TITLE).txt data/uma_pe_equipartition_ratio_vs_dt_$(TITLE).txt data/config_pe_equipartition_ratio_vs_dt_$(TITLE).txt

########################### THESIS STUFF #################################

THESIS_FIGURES = plots/trajectory-plot_thesis.pdf plots/stepping_length_histogram_thesis.pdf plots/stepping_time_histogram_thesis.pdf

thesis_stuff.pdf: thesis_stuff/thesis_stuff.tex thesis_stuff/thesis_stuff.bib $(FIGURES)
	cd thesis_stuff && xelatex -interaction nonstopmode -halt-on-error thesis_stuff.tex && bibtex thesis_stuff && xelatex -interaction nonstopmode -halt-on-error thesis_stuff.tex && xelatex -interaction nonstopmode -halt-on-error thesis_stuff.tex && mv thesis_stuff.pdf ..

thesis: $(THESIS_FIGURES)
	cd thesis && $(MAKE)

clean:
	rm -f *.txt
	rm -f *.o
	rm -f plot
	rm -f test_onebound
	rm -f test_bothbound
	rm -f ob_PE_correlation_vs_time
	rm -f ob_PE_equipartition_ratio_vs_time
	rm -f ob_PE_equipartition_ratio_average_vs_time
	rm -f ob_PE_equipartition_ratio_average_vs_force_ratio
	rm -f ob_PE_equipartition_ratio_average_vs_spring_constant
	rm -f ob_PE_equipartition_ratio_average_vs_temperature
	rm -f ob_PE_equipartition_ratio_average_vs_dt
	rm -f ob_PE_log_equipartition_vs_log_iterations
	rm -f onebound_PE_equipartition_correlation
	rm -f ob_equipartition_test
	rm -f generate_stepping_data
	rm -f data.txt
	rm -f config.txt
	rm -f latex/*.aux
	rm -f latex/*.log
	rm -f latex/*.pdf
	rm -f latex/latexlog.txt
	rm -f *~
	rm -f simulations/*~
	rm -f *#
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm data/thesis_data.txt
	rm data/thesis_plot.txt
	cd figures && $(MAKE) clean
	cd thesis && $(MAKE) clean
