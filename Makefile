CPPFLAGS = -ggdb

all: walk data.txt

Motion_Equations/Leftbound_Solution.txt: Mathematica/Leftbound_Solution.m
	math -noprompt -script Mathematica/Leftbound_Solution.m
	
Motion_Equations/Rightbound_Solution.txt: Mathematica/Rightbound_Solution.m
	math -noprompt -script Mathematica/Rightbound_Solution.m
	
Motion_Equations/Bothbound_Solution.txt: Mathematica/Bothbound_Solution.m
	math -noprompt -script Mathematica/Bothbound_Solution.m

dynein_motion_functions.cpp: replace.py Motion_Equations/Leftbound_Solution.txt 
	python replace.py

walk: dynein_walk.cpp dynein_struct.cpp utilities.cpp dynein_struct.h dynein_motion_functions.cpp
	g++ dynein_walk.cpp dynein_struct.cpp dynein_motion_functions.cpp utilities.cpp -o walk $(CPPFLAGS)

test: dynein_ftest.cpp dynein_struct.cpp dynein_motion_functions.cpp dynein_struct.h dynein_motion_functions.cpp
	g++ dynein_ftest.cpp dynein_struct.cpp dynein_motion_functions.cpp utilities.cpp -o ftest $(CPPFLAGS)

clean:
	rm -f *.o
	rm -f walk
	rm -f ftest
	rm -f data.txt
	rm -f config.txt 
