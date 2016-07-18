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

bothbound_equipartition_test: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o  simulations/bothbound_equipartition_test.cpp FORCE
	g++ -c simulations/bothbound_equipartition_test.cpp $(CPPFLAGS)
	g++ bothbound_equipartition_test.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o bothbound_equipartition_test
	./bothbound_equipartition_test

ob_equipartition_test: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_equipartition_test.cpp FORCE
	g++ -c simulations/ob_equipartition_test.cpp $(CPPFLAGS)
	g++ ob_equipartition_test.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_equipartition_test
	./ob_equipartition_test

ob_PE_correlation_vs_time: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_correlation_vs_time.cpp FORCE
	g++ -c simulations/ob_PE_correlation_vs_time.cpp $(CPPFLAGS)
	g++ ob_PE_correlation_vs_time.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_correlation_vs_time

ob_PE_correlation_vs_time_plot: ob_PE_correlation_vs_time FORCE
	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
	mkdir -p plots
	mkdir -p data
	./ob_PE_correlation_vs_time $(TITLE)
	./make_plot.py --figtitle="$(TITLE)" --xlabel="Tau (s)" --ylabel="Correlation" data/bba_pe_correlation_$(TITLE).txt data/bma_pe_correlation_$(TITLE).txt data/ta_pe_correlation_$(TITLE).txt data/uma_pe_correlation_$(TITLE).txt

ob_PE_equipartition_ratio_average_vs_time: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_equipartition_ratio_average_vs_time.cpp FORCE
	g++ -c simulations/ob_PE_equipartition_ratio_average_vs_time.cpp $(CPPFLAGS)
	g++ ob_PE_equipartition_ratio_average_vs_time.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_average_vs_time

ob_PE_equipartition_ratio_average_vs_time_plot: ob_PE_equipartition_ratio_average_vs_time FORCE
	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
	mkdir -p plots
	mkdir -p data
	./ob_PE_equipartition_ratio_average_vs_time $(TITLE)
	./make_plot.py --figtitle="$(TITLE)" --xlabel="Runtime (s)" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/bba_pe_equipartition_ratio_v_time_$(TITLE).txt data/bma_pe_equipartition_ratio_v_time_$(TITLE).txt data/ta_pe_equipartition_ratio_v_time_$(TITLE).txt data/uma_pe_equipartition_ratio_v_time_$(TITLE).txt

ob_PE_equipartition_ratio_vs_time: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_equipartition_ratio_vs_time.cpp FORCE
	g++ -c simulations/ob_PE_equipartition_ratio_vs_time.cpp $(CPPFLAGS)
	g++ ob_PE_equipartition_ratio_vs_time.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_vs_time

ob_PE_equipartition_ratio_vs_time_plot: ob_PE_equipartition_ratio_vs_time FORCE
	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
	mkdir -p plots
	mkdir -p data
	./ob_PE_equipartition_ratio_vs_time $(TITLE)
	./make_plot.py --figtitle="$(TITLE)" --xlabel="Runtime (s)" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/bba_pe_equipartition_ratio_$(TITLE).txt data/bma_pe_equipartition_ratio_$(TITLE).txt data/ta_pe_equipartition_ratio_$(TITLE).txt data/uma_pe_equipartition_ratio_$(TITLE).txt

ob_PE_equipartition_ratio_average_vs_spring_constant: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_equipartition_ratio_average_vs_spring_constant.cpp FORCE
	g++ -c simulations/ob_PE_equipartition_ratio_average_vs_spring_constant.cpp $(CPPFLAGS)
	g++ ob_PE_equipartition_ratio_average_vs_spring_constant.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_average_vs_spring_constant

ob_PE_equipartition_ratio_average_vs_spring_constant_plot: ob_PE_equipartition_ratio_average_vs_spring_constant FORCE
	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
	mkdir -p plots
	mkdir -p data
	./ob_PE_equipartition_ratio_average_vs_spring_constant $(TITLE)
	./make_plot.py --figtitle="$(TITLE)" --logx --xlabel="Spring constant (nm^2*kg/s^2)" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/bba_pe_equipartition_ratio_vs_c_$(TITLE).txt data/bma_pe_equipartition_ratio_vs_c_$(TITLE).txt data/ta_pe_equipartition_ratio_vs_c_$(TITLE).txt data/uma_pe_equipartition_ratio_vs_c_$(TITLE).txt

ob_PE_equipartition_ratio_average_vs_force_ratio: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_equipartition_ratio_average_vs_force_ratio.cpp FORCE
	g++ -c simulations/ob_PE_equipartition_ratio_average_vs_force_ratio.cpp $(CPPFLAGS)
	g++ ob_PE_equipartition_ratio_average_vs_force_ratio.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_average_vs_force_ratio

ob_PE_equipartition_ratio_average_vs_force_ratio_plot: ob_PE_equipartition_ratio_average_vs_force_ratio FORCE
	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
	mkdir -p plots
	mkdir -p data
	./ob_PE_equipartition_ratio_average_vs_force_ratio $(TITLE)
	./make_plot.py --scatter --figtitle="$(TITLE)" --xlabel="Brownian / conformational force variance" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/bba_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt data/bma_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt data/ta_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt data/uma_pe_equipartition_ratio_vs_f_ratio_$(TITLE).txt

ob_PE_equipartition_ratio_average_vs_temperature: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_equipartition_ratio_average_vs_temperature.cpp FORCE
	g++ -c simulations/ob_PE_equipartition_ratio_average_vs_temperature.cpp $(CPPFLAGS)
	g++ ob_PE_equipartition_ratio_average_vs_temperature.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_average_vs_temperature

ob_PE_equipartition_ratio_average_vs_temperature_plot: ob_PE_equipartition_ratio_average_vs_temperature FORCE
	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
	mkdir -p plots
	mkdir -p data
	./ob_PE_equipartition_ratio_average_vs_temperature $(TITLE)
	./make_plot.py --figtitle="$(TITLE)" --xlabel="Temp (K)" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/bba_pe_equipartition_ratio_vs_T_$(TITLE).txt data/bma_pe_equipartition_ratio_vs_T_$(TITLE).txt data/ta_pe_equipartition_ratio_vs_T_$(TITLE).txt data/uma_pe_equipartition_ratio_vs_T_$(TITLE).txt

ob_PE_log_equipartition_vs_log_iterations: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_log_equipartition_vs_log_iterations.cpp FORCE
	g++ -c simulations/ob_PE_log_equipartition_vs_log_iterations.cpp $(CPPFLAGS)
	g++ ob_PE_log_equipartition_vs_log_iterations.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_log_equipartition_vs_log_iterations $(LIBRARIES)

ob_PE_log_equipartition_vs_log_iterations_plot: ob_PE_log_equipartition_vs_log_iterations FORCE
	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
	mkdir -p plots
	mkdir -p data
	./ob_PE_log_equipartition_vs_log_iterations "$(TITLE)"
	./make_plot.py --loglog --figtitle="$(TITLE)" --xlabel="log(iterations)" --ylabel="log(| PE / ET | - 1)" --hline=1.0 data/bba_pe_log_equipartition_vs_log_iterations_$(TITLE).txt data/bma_pe_log_equipartition_vs_log_iterations_$(TITLE).txt data/ta_pe_log_equipartition_vs_log_iterations_$(TITLE).txt data/uma_pe_log_equipartition_vs_log_iterations_$(TITLE).txt

ob_PE_equipartition_ratio_average_vs_dt: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/ob_PE_equipartition_ratio_average_vs_dt.cpp FORCE
	g++ -c simulations/ob_PE_equipartition_ratio_average_vs_dt.cpp $(CPPFLAGS)
	g++ ob_PE_equipartition_ratio_average_vs_dt.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o ob_PE_equipartition_ratio_average_vs_dt

ob_PE_equipartition_ratio_average_vs_dt_plot: ob_PE_equipartition_ratio_average_vs_dt FORCE
	@echo "\nUse TITLE='yourtitle' to give plot a title\n"
	mkdir -p plots
	mkdir -p data
	./ob_PE_equipartition_ratio_average_vs_dt $(TITLE)
	./make_plot.py --figtitle="$(TITLE)" --xlabel="log(dt)" --ylabel="PE / 0.5*kb*T" --logx data/bba_pe_equipartition_ratio_vs_dt_$(TITLE).txt data/bma_pe_equipartition_ratio_vs_dt_$(TITLE).txt data/ta_pe_equipartition_ratio_vs_dt_$(TITLE).txt data/uma_pe_equipartition_ratio_vs_dt_$(TITLE).txt

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
