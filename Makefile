CPPFLAGS = -ggdb

all: walk data.txt

solve: Mathematica/*
	math -noprompt -script Mathematica/Dyn_br_left_solve.m
	math -noprompt -script Mathematica/Dyn_br_both_solve.m
	math -noprompt -script Mathematica/Dyn_br_right_solve.m

dynein_motion_functions.cpp: replace.py Motion_Equations/DyneinBrownianBothboundSolutionsUnsimplified.txt \
				Motion_Equations/DyneinBrownianLeftboundSolutionsUnsimplified.txt \
				Motion_Equations/DyneinBrownianRightboundSolutionsUnsimplified.txt
	python replace.py

walk: dynein_walk.cpp dynein_struct.cpp utilities.cpp dynein_struct.h dynein_motion_functions.cpp
	g++ dynein_walk.cpp dynein_struct.cpp dynein_motion_functions.cpp utilities.cpp -o walk $(CPPFLAGS)

data.txt: walk
	time ./walk

test: dynein_ftest.cpp dynein_struct.cpp dynein_motion_functions.cpp dynein_struct.h dynein_motion_functions.cpp
	g++ dynein_ftest.cpp dynein_struct.cpp dynein_motion_functions.cpp utilities.cpp -o ftest $(CPPFLAGS)
	math -noprompt -script Mathematica/Dyn_br_left_test.m
	math -noprompt -script Mathematica/Dyn_br_both_test.m
	math -noprompt -script Mathematica/Dyn_br_right_test.m

#math -noprompt -script Mathematica/Dyn_br_both_test.m

clean:
	rm -f *.o
	rm -f walk
	rm -f ftest
	rm -f data.txt
	rm -f config.txt

#	rm -f prevents rm from complaining when there is no such file to delete
