CPPFLAGS = -std=c++11 -g -Wall -Werror -O2

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

bothbound_equipartition_test: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o  simulations/bothbound_equipartition_test.cpp FORCE
	g++ -c simulations/bothbound_equipartition_test.cpp $(CPPFLAGS)
	g++ bothbound_equipartition_test.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o bothbound_equipartition_test
	./bothbound_equipartition_test

onebound_equipartition_test: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/onebound_equipartition_test.cpp FORCE
	g++ -c simulations/onebound_equipartition_test.cpp $(CPPFLAGS)
	g++ onebound_equipartition_test.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o onebound_equipartition_test
	./onebound_equipartition_test

onebound_PE_equipartition_correlation_plot: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o simulations/onebound_PE_equipartition_correlation.cpp simulations/simulations.h FORCE
	mkdir -p movies
	g++ -c simulations/onebound_PE_equipartition_correlation.cpp $(CPPFLAGS)
	g++ onebound_PE_equipartition_correlation.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations.o -o onebound_PE_equipartition_correlation_plot
	#./onebound_PE_equipartition_correlation_plot
	#./make_plot.py --figtitle="Equipartition ratio for conformational PEs" --xlabel="Runtime (s)" --ylabel="PE / 0.5*kb*T" --hline=1.0 bba_pe_equipartition_ratio.txt bma_pe_equipartition_ratio.txt ta_pe_equipartition_ratio.txt uma_pe_equipartition_ratio.txt
	#./make_plot.py --figtitle="Correlation function for conformational PEs" --xlabel="Tau (s)" --ylabel="Correlation" bba_pe_correlation.txt bma_pe_correlation.txt ta_pe_correlation.txt uma_pe_correlation.txt

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

thesis_stuff/thesis_stuff.pdf: thesis_stuff/thesis_stuff.tex
	cd thesis_stuff && xelatex thesis_stuff.tex

thesis_stuff/log.pdf: thesis_stuff/random_derivations.tex
	cd random_derivations && xelatex random_derivations.tex

clean:
	rm -f *.txt
	rm -f *.o
	rm -f plot
	rm -f test_onebound
	rm -f test_bothbound
	rm -f PE_correlation_function
	rm -f data.txt
	rm -f config.txt
	rm -f latex/*.aux
	rm -f latex/*.log
	rm -f latex/*.pdf
	rm -f latex/latexlog.txt
	rm -f *~
	rm -f *#

FORCE:
