CPPFLAGS = -ggdb -Wall

all: walk data.txt

walk: dynein_walk.cpp dynein_struct.cpp utilities.cpp dynein_struct.h
	lint dynein_walk.cpp dynein_struct.cpp utilities.cpp dynein_struct.h
	g++ dynein_walk.cpp dynein_struct.cpp utilities.cpp -o walk $(CPPFLAGS)

test: dynein_ftest.cpp dynein_struct.cpp utilities.cpp dynein_struct.h
	g++ dynein_ftest.cpp dynein_struct.cpp utilities.cpp -o ftest $(CPPFLAGS)

clean:
	rm -f *.o
	rm -f walk
	rm -f ftest
	rm -f data.txt
	rm -f config.txt
