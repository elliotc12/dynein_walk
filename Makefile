CPPFLAGS = -g -Wall -O2

all: test walk data.txt derivation.pdf

dynein_walk.o: dynein_walk.cpp dynein_struct.h
	g++ -c dynein_walk.cpp $(CPPFLAGS)

dynein_struct.o: dynein_struct.cpp dynein_struct.h
	g++ -c dynein_struct.cpp $(CPPFLAGS)

dynein_test.o: dynein_test.cpp dynein_struct.h
	g++ -c dynein_test.cpp $(CPPFLAGS)

utilities.o: utilities.cpp dynein_struct.h
	g++ -c utilities.cpp $(CPPFLAGS)

walk: test dynein_walk.o dynein_struct.o utilities.o
	./test
	g++ dynein_walk.o dynein_struct.o utilities.o -o walk

test: dynein_test.o dynein_struct.o utilities.o
	g++ dynein_test.o dynein_struct.o utilities.o -o test

data.txt: walk
	./walk

plot: data.txt
	./plot.py

derivation.pdf: latex/derivation.tex
	cd latex && pdflatex derivation.tex && mv derivation.pdf ..

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
