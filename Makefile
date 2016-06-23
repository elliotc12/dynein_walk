CPPFLAGS = -std=c++11 -g -Wall -Werror -O2

all: plot

derivation.pdf: latex/derivation.tex
	cd latex && pdflatex derivation.tex && mv derivation.pdf ..

derivation_confirmation.pdf: latex/derivation_confirmation.tex
	cd latex && pdflatex derivation_confirmation.tex && mv derivation_confirmation.pdf ..

test_bothbound.o: test_bothbound.cpp dynein_struct.h
	g++ -c test_bothbound.cpp $(CPPFLAGS)

dynein_struct_onebound.o: dynein_struct_onebound.cpp dynein_struct.h
	g++ -c dynein_struct_onebound.cpp $(CPPFLAGS)

dynein_struct_bothbound.o: dynein_struct_bothbound.cpp dynein_struct.h
	g++ -c dynein_struct_bothbound.cpp $(CPPFLAGS)

dynein_simulate.o: dynein_simulate.cpp dynein_struct_onebound.cpp dynein_struct_bothbound.cpp
	g++ -c dynein_simulate.cpp $(CPPFLAGS)

equipartition_test: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o test_onebound test_bothbound simulations/equipartition_test.cpp
	g++ -c simulations/equipartition_test.cpp $(CPPFLAGS)
	g++ equipartition_test.o dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o -o equipartition_test
	./equipartition_test

test_onebound.o: test_onebound.cpp dynein_struct.h
	g++ -c test_onebound.cpp $(CPPFLAGS)

figures: figures/*
	cd figures && make

plot: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations/plot.cpp test_onebound test_bothbound FORCE
	g++ -c simulations/plot.cpp $(CPPFLAGS) -o plot.o
	g++ dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o plot.o -o plot
	./simulate.py veryfast verylong natural $(OPT)

paper.pdf: latex/paper.tex
	cd latex && pdflatex paper.tex && mv paper.pdf ..

save: dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o simulations/plot.cpp test_onebound test_bothbound
	mkdir -p movies
	mkdir -p PNGs
	g++ -c simulations/plot.cpp $(CPPFLAGS) -o plot.o
	g++ dynein_simulate.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o plot.o -o plot
	./simulate.py veryfast verylong natural save $(NAME)

test_bothbound: test_bothbound.o dynein_struct_bothbound.o dynein_struct_onebound.o utilities.o
	g++ test_bothbound.o dynein_struct_bothbound.o dynein_struct_onebound.o utilities.o -o test_bothbound
	./test_bothbound

test_onebound: test_onebound.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o
	g++ test_onebound.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o -o test_onebound
	./test_onebound

utilities.o: utilities.cpp dynein_struct.h
	g++ -c utilities.cpp $(CPPFLAGS)

clean:
	rm -f *.o
	rm -f plot
	rm -f test_onebound
	rm -f test_bothbound
	rm -f data.txt
	rm -f config.txt
	rm -f latex/*.aux
	rm -f latex/*.log
	rm -f latex/*.pdf
	rm -f latex/latexlog.txt
	rm -f *~
	rm -f *#

FORCE:
