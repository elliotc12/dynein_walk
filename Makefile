CPPFLAGS = -ggdb

solve: Mathematica/*
	math -noprompt -script Mathematica/Dyn_br_left_solve.m
	math -noprompt -script Mathematica/Dyn_br_both_solve.m
	math -noprompt -script Mathematica/Dyn_br_right_solve.m

walk: replace.py Motion_Equations/* dynein_walk.cpp dynein_struct.cpp utilities.cpp dynein_struct.h
	python replace.py
	g++ dynein_walk.cpp dynein_struct.cpp dynein_motion_functions.cpp utilities.cpp -o walk $(CPPFLAGS)
	
test: replace.py Motion_Equations/* dynein_ftest.cpp dynein_struct.cpp dynein_motion_functions.cpp dynein_struct.h
	python replace.py
	g++ dynein_ftest.cpp dynein_struct.cpp dynein_motion_functions.cpp -o ftest $(CPPFLAGS)

clean:
	rm -f *.o
	rm -f walk
	rm -f ftest
	rm -f data.txt
	rm -f config.txt

#	rm -f prevents rm from complaining when there is no such file to delete
