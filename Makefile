CPPFLAGS = -g -Wall -O2

all: walk data.txt

dynein_walk.o: dynein_walk.cpp dynein_struct.h
	#lint dynein_walk.cpp
	g++ -c dynein_walk.cpp $(CPPFLAGS)

dynein_struct.o: dynein_struct.cpp dynein_struct.h
	#lint dynein_struct.cpp
	g++ -c dynein_struct.cpp $(CPPFLAGS)

dynein_ftest.o: dynein_ftest.cpp dynein_struct.h
	#lint dynein_ftest.cpp
	g++ -c dynein_ftest.cpp $(CPPFLAGS)

utilities.o: utilities.cpp dynein_struct.h
	#lint utilities.cpp
	g++ -c utilities.cpp $(CPPFLAGS)

walk: dynein_walk.o dynein_struct.o utilities.o
	g++ dynein_walk.o dynein_struct.o utilities.o -o walk

test: dynein_ftest.o dynein_struct.o utilities.o
	g++ dynein_ftest.o dynein_struct.o utilities.o -o ftest

data.txt: walk
	./walk

clean:
	rm -f *.o
	rm -f walk
	rm -f ftest
	rm -f data.txt
	rm -f config.txt
