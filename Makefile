CPPFLAGS = -std=c++11 -g -Wall -Werror -O2

#all: test walk plot
all: walk plot

dynein_walk.o: dynein_walk.cpp dynein_struct.h
	g++ -c dynein_walk.cpp $(CPPFLAGS)

dynein_struct_onebound.o: dynein_struct_onebound.cpp dynein_struct.h
	g++ -c dynein_struct_onebound.cpp $(CPPFLAGS)

dynein_struct_bothbound.o: dynein_struct_bothbound.cpp dynein_struct.h
	g++ -c dynein_struct_bothbound.cpp $(CPPFLAGS)

dynein_test.o: dynein_test.cpp dynein_struct.h
	g++ -c dynein_test.cpp $(CPPFLAGS)

utilities.o: utilities.cpp dynein_struct.h
	g++ -c utilities.cpp $(CPPFLAGS)

#walk: test dynein_walk.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o
walk: dynein_walk.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o
#	./test
	g++ dynein_walk.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o -o walk

test: dynein_test.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o
	g++ dynein_test.o dynein_struct_onebound.o dynein_struct_bothbound.o utilities.o -o test

plot: walk
	./simulate.py veryfast short natural loop

derivation.pdf: latex/derivation.tex
	cd latex && pdflatex derivation.tex && mv derivation.pdf ..

derivation_confirmation.pdf: latex/derivation_confirmation.tex
	cd latex && pdflatex derivation_confirmation.tex && mv derivation_confirmation.pdf ..

paper.pdf: latex/paper.tex
	cd latex && pdflatex paper.tex && mv paper.pdf ..

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
