CPPFLAGS = -std=c++11 -g -Wall -Werror -O2

all: walk plot

derivation.pdf: latex/derivation.tex
	cd latex && pdflatex derivation.tex && mv derivation.pdf ..

derivation_confirmation.pdf: latex/derivation_confirmation.tex
	cd latex && pdflatex derivation_confirmation.tex && mv derivation_confirmation.pdf ..

dynein_walk.o: dynein_walk.cpp dynein_struct.h
	g++ -c dynein_walk.cpp $(CPPFLAGS)

test_bothbound.o: test_bothbound.cpp dynein_struct.h
	g++ -c test_bothbound.cpp $(CPPFLAGS)

dynein_struct_onebound.o: dynein_struct_onebound.cpp dynein_struct.h
	g++ -c dynein_struct_onebound.cpp $(CPPFLAGS)

dynein_struct_bothbound.o: dynein_struct_bothbound.cpp dynein_struct.h
	g++ -c dynein_struct_bothbound.cpp $(CPPFLAGS)

test_onebound.o: test_onebound.cpp dynein_struct.h
	g++ -c test_onebound.cpp $(CPPFLAGS)

figures: figures/*
	cd figures && make

plot: walk
	./simulate.py veryfast verylong natural

plot-step: walk
	./simulate.py veryfast normal natural step

paper.pdf: latex/paper.tex
	cd latex && pdflatex paper.tex && mv paper.pdf ..

save: walk
	mkdir -p movies
	mkdir -p PNGs
	./simulate.py veryfast normal natural save $(NAME)

test_bothbound: test_bothbound.o dynein_struct_bothbound.o dynein_struct_onebound.o utilities.o
	g++ test_bothbound.o dynein_struct_bothbound.o dynein_struct_onebound.o utilities.o -o test_bothbound
	./test_bothbound

test_onebound: test_onebound.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o
	g++ test_onebound.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o -o test_onebound
	./test_onebound

walk: dynein_walk.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o test_onebound test_bothbound
	g++ dynein_walk.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o -o walk

utilities.o: utilities.cpp dynein_struct.h
	g++ -c utilities.cpp $(CPPFLAGS)

clean:
	rm -f *.o
	rm -f walk
	rm -f test
	rm -f data.txt
	rm -f config.txt
	rm -f latex/*.aux
	rm -f latex/*.log
	rm -f latex/*.pdf
	rm -f latex/latexlog.txt

FORCE:
