CPPFLAGS = -std=c++11 -g -Wall -Werror -O2
LIBRARIES = -lm

.PHONY: test_bothbound test_onebound clean plot

all: test_onebound test_bothbound plot

derivation.pdf: latex/derivation.tex
	cd latex && pdflatex derivation.tex && mv derivation.pdf ..

derivation_confirmation.pdf: latex/derivation_confirmation.tex
	cd latex && pdflatex derivation_confirmation.tex && mv derivation_confirmation.pdf ..

test_bothbound.o: test_bothbound.cpp dynein_struct.h default_parameters.h
	g++ -c test_bothbound.cpp $(CPPFLAGS)

dynein_struct_onebound.o: dynein_struct_onebound.cpp dynein_struct.h default_parameters.h
	g++ -c dynein_struct_onebound.cpp $(CPPFLAGS)

dynein_struct_bothbound.o: dynein_struct_bothbound.cpp dynein_struct.h default_parameters.h
	g++ -c dynein_struct_bothbound.cpp $(CPPFLAGS)

dynein_simulate.o: dynein_simulate.cpp dynein_struct_onebound.cpp dynein_struct_bothbound.cpp default_parameters.h
	g++ -c dynein_simulate.cpp $(CPPFLAGS)

simulations.o: simulations/simulations.cpp dynein_struct.h default_parameters.h
	g++ -c simulations/simulations.cpp $(CPPFLAGS) -o simulations.o

test_onebound.o: test_onebound.cpp dynein_struct.h default_parameters.h
	g++ -c test_onebound.cpp $(CPPFLAGS)

figures: figures/*
	cd figures && make

plot: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations/plot.cpp
	@echo "CONF = natural / pretty"
	g++ -c simulations/plot.cpp $(CPPFLAGS) -o plot.o
	g++ dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o plot.o -o plot
	./simulate.py veryfast verylong $(CONF) $(OPT)

paper.pdf: latex/paper.tex
	cd latex && pdflatex paper.tex && mv paper.pdf ..

save: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations/plot.cpp
	mkdir -p movies
	mkdir -p PNGs
	g++ -c simulations/plot.cpp $(CPPFLAGS) -o plot.o
	g++ dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o plot.o -o plot
	./simulate.py veryfast verylong natural save $(NAME)

test_bothbound: test_bothbound.o dynein_struct_bothbound.o dynein_struct_onebound.o utilities.o dynein_simulate.o
	g++ test_bothbound.o dynein_struct_bothbound.o dynein_struct_onebound.o dynein_simulate.o utilities.o -o test_bothbound
	./test_bothbound

test_onebound: test_onebound.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o dynein_simulate.o
	g++ test_onebound.o dynein_struct_onebound.o dynein_struct_bothbound.o dynein_simulate.o utilities.o -o test_onebound
	./test_onebound

utilities.o: utilities.cpp dynein_struct.h default_parameters.h
	g++ -c utilities.cpp $(CPPFLAGS)

######################### SIMULATION STUFF ###############################
TITLE = defaultplot

create_ob_plots: simulations/create_ob_plots.cpp dynein_struct.h default_parameters.h simulations/simulation_defaults.h dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o
	g++ -c simulations/create_ob_plots.cpp $(CPPFLAGS)
	g++ dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o create_ob_plots.o -o create_ob_plots

ob_plots: create_ob_plots FORCE
	@echo "Use TITLE='yourtitle' to give plot a title"
	./create_ob_plots $(TITLE)
	mkdir -p plots
	./make_plot.py --figtitle="Correlation_function_$(TITLE)" --xlabel="Tau (s)" --ylabel="Correlation" data/ob_bba_pe_$(TITLE)_correlation_fn.txt data/ob_bma_pe_$(TITLE)_correlation_fn.txt data/ob_ta_pe_$(TITLE)_correlation_fn.txt data/ob_uma_pe_$(TITLE)_correlation_fn.txt data/ob_config_$(TITLE).txt
	./make_plot.py --figtitle="Locally averaged PE_vs_time_$(TITLE)" --skiprows=100 --xlabel="Runtime (s)" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/ob_bba_pe_$(TITLE).txt data/ob_bma_pe_$(TITLE).txt data/ob_ta_pe_$(TITLE).txt data/ob_uma_pe_$(TITLE).txt data/ob_config_$(TITLE).txt
	./make_plot.py --figtitle="PE_average_vs_time_$(TITLE)" --xlabel="Runtime (s)" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/ob_bba_pe_$(TITLE)_eq_ave.txt data/ob_bma_pe_$(TITLE)_eq_ave.txt data/ob_ta_pe_$(TITLE)_eq_ave.txt data/ob_uma_pe_$(TITLE)_eq_ave.txt data/ob_config_$(TITLE).txt
	./make_plot.py --logx --logy --figtitle="Log_error_vs_log_time_$(TITLE)" --xlabel="log(iterations)" --ylabel="log(| PE / ET - 1|)" --hline=1.0 data/ob_bba_pe_$(TITLE)_log_error.txt data/ob_bma_pe_$(TITLE)_log_error.txt data/ob_ta_pe_$(TITLE)_log_error.txt data/ob_uma_pe_$(TITLE)_log_error.txt data/ob_config_$(TITLE).txt
	./make_plot.py --figtitle="Angle_vs_time_$(TITLE)" --xlabel="Runtime (s)" --ylabel="Angle" data/ob_bba_angle_$(TITLE).txt data/ob_bma_angle_$(TITLE).txt data/ob_ta_angle_$(TITLE).txt data/ob_uma_angle_$(TITLE).txt data/ob_config_$(TITLE).txt

generate_onebound_data: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/generate_onebound_data.cpp default_parameters.h dynein_struct.h FORCE
	mkdir -p data
	g++ -c simulations/generate_onebound_data.cpp $(CPPFLAGS)
	g++ generate_onebound_data.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o generate_onebound_data
	./generate_onebound_data $(TITLE)

bothbound_equipartition_test: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o  simulations/bothbound_equipartition_test.cpp FORCE
	g++ -c simulations/bothbound_equipartition_test.cpp $(CPPFLAGS)
	g++ bothbound_equipartition_test.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o bothbound_equipartition_test
	./bothbound_equipartition_test

ob_equipartition_test: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_equipartition_test.cpp FORCE
	g++ -c simulations/ob_equipartition_test.cpp $(CPPFLAGS)
	g++ ob_equipartition_test.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_equipartition_test
	./ob_equipartition_test

ob_PE_equipartition_ratio_average_vs_spring_constant: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_equipartition_ratio_average_vs_spring_constant.cpp FORCE
	g++ -c simulations/ob_PE_equipartition_ratio_average_vs_spring_constant.cpp $(CPPFLAGS)
	g++ ob_PE_equipartition_ratio_average_vs_spring_constant.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_average_vs_spring_constant

ob_PE_equipartition_ratio_average_vs_spring_constant_plot: ob_PE_equipartition_ratio_average_vs_spring_constant FORCE
	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
	mkdir -p plots
	mkdir -p data
	./ob_PE_equipartition_ratio_average_vs_spring_constant $(TITLE)
	./make_plot.py --figtitle="$(TITLE)" --logx --xlabel="Spring constant (nm^2*kg/s^2)" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/bba_pe_equipartition_ratio_vs_c_$(TITLE).txt data/bma_pe_equipartition_ratio_vs_c_$(TITLE).txt data/ta_pe_equipartition_ratio_vs_c_$(TITLE).txt data/uma_pe_equipartition_ratio_vs_c_$(TITLE).txt data/config_pe_equipartition_ratio_vs_c_$(TITLE).txt

ob_PE_equipartition_ratio_average_vs_force_ratio: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_equipartition_ratio_average_vs_force_ratio.cpp FORCE
	g++ -c simulations/ob_PE_equipartition_ratio_average_vs_force_ratio.cpp $(CPPFLAGS)
	g++ ob_PE_equipartition_ratio_average_vs_force_ratio.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_average_vs_force_ratio

ob_PE_equipartition_ratio_average_vs_force_ratio_plot: ob_PE_equipartition_ratio_average_vs_force_ratio FORCE
	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
	mkdir -p plots
	mkdir -p data
	./ob_PE_equipartition_ratio_average_vs_force_ratio $(TITLE)
	./make_plot.py --scatter --figtitle="$(TITLE)" --xlabel="Brownian / conformational force variance" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/bba_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt data/bma_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt data/ta_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt data/uma_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt data/config_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt

ob_PE_equipartition_ratio_average_vs_temperature: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_equipartition_ratio_average_vs_temperature.cpp FORCE
	g++ -c simulations/ob_PE_equipartition_ratio_average_vs_temperature.cpp $(CPPFLAGS)
	g++ ob_PE_equipartition_ratio_average_vs_temperature.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_average_vs_temperature

ob_PE_equipartition_ratio_average_vs_temperature_plot: ob_PE_equipartition_ratio_average_vs_temperature FORCE
	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
	mkdir -p plots
	mkdir -p data
	./ob_PE_equipartition_ratio_average_vs_temperature $(TITLE)
	./make_plot.py --figtitle="$(TITLE)" --xlabel="Temp (K)" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/bba_pe_equipartition_ratio_vs_T_$(TITLE).txt data/bma_pe_equipartition_ratio_vs_T_$(TITLE).txt data/ta_pe_equipartition_ratio_vs_T_$(TITLE).txt data/uma_pe_equipartition_ratio_vs_T_$(TITLE).txt data/config_pe_equipartition_ratio_vs_T_$(TITLE).txt

ob_PE_equipartition_ratio_average_vs_dt: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_equipartition_ratio_average_vs_dt.cpp FORCE
	g++ -c simulations/ob_PE_equipartition_ratio_average_vs_dt.cpp $(CPPFLAGS)
	g++ ob_PE_equipartition_ratio_average_vs_dt.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_average_vs_dt

ob_PE_equipartition_ratio_average_vs_dt_plot: ob_PE_equipartition_ratio_average_vs_dt FORCE
	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
	mkdir -p plots
	mkdir -p data
	./ob_PE_equipartition_ratio_average_vs_dt $(TITLE)
	./make_plot.py --figtitle="$(TITLE)" --xlabel="log(dt)" --ylabel="PE / 0.5*kb*T" --logx data/bba_pe_equipartition_ratio_vs_dt_$(TITLE).txt data/bma_pe_equipartition_ratio_vs_dt_$(TITLE).txt data/ta_pe_equipartition_ratio_vs_dt_$(TITLE).txt data/uma_pe_equipartition_ratio_vs_dt_$(TITLE).txt data/config_pe_equipartition_ratio_vs_dt_$(TITLE).txt

########################### THESIS STUFF #################################

thesis_stuff/thesis_stuff.pdf: thesis_stuff/thesis_stuff.tex FORCE
	cd thesis_stuff && xelatex thesis_stuff.tex

thesis_stuff/random_derivations.pdf: thesis_stuff/random_derivations.tex FORCE
	cd thesis_stuff && xelatex random_derivations.tex

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
	rm -f data.txt
	rm -f config.txt
	rm -f latex/*.aux
	rm -f latex/*.log
	rm -f latex/*.pdf
	rm -f latex/latexlog.txt
	rm -f *~
	rm -f *#

FORCE:
