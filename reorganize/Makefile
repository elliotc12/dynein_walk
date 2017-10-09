CPPFLAGS = -std=c++98 -g -Werror -O2 -Wall
LIBRARIES = -lm

.PHONY: clean

.PRECIOUS: data/stepping_data_%.txt data/stepping_config_%.txt data/stepping_movie_data_%.txt data/bothbound_data_%.bin data/onebound_data_%.bin data/ob_config_%.txt data/bb_config_%.txt # prevent nonexistant data files from being deleted after creation+use

PAPERS = papers/fft_speedup/convolution_theorem.pdf
DRAW = scripts/dynein/draw/motor_domain.py scripts/dynein/draw/tail.py

all: generate_stepping_data $(PAPERS) $(DRAW)

src/version-info.h: SHELL:=/bin/bash
src/version-info.h: ../.git/refs/heads/master $(wildcard src/*.cpp) Makefile
	UNAMESTR=$$(uname); if [[ "$$UNAMESTR" == 'Linux' ]]; then \
	echo -n static const char '*version' = '"' > src/version-info.h; \
	git describe --dirty --tags | tr -d '\n' >> src/version-info.h; \
	echo -n '-=-' >> src/version-info.h; \
	date -Ins | tr -d '\n' >> src/version-info.h; \
	echo '";' >> src/version-info.h; \
	elif [[ "$$UNAMESTR" == 'Darwin' ]]; then \
	    echo 'static const char *version = "mac-breaks-version-stuff";' > src/version-info.h; \
	fi;

src/generate_stepping_data.o: src/generate_stepping_data.cpp src/version-info.h

generate_stepping_data: src/generate_stepping_data.o src/dynein_simulate.o \
			src/dynein_struct_onebound.o src/dynein_struct_bothbound.o src/utilities.o
	$(CXX) -o generate_stepping_data $^

papers/fft_speedup/convolution_theorem.pdf: papers/fft_speedup/convolution_theorem.tex
	cd papers/fft_speedup && pdflatex convolution_theorem &&  pdflatex convolution_theorem

## draw module stuff ##
scripts/dynein/draw/motor_domain.py: scripts/dynein/draw/create_MD_array.py scripts/dynein/draw/outer_coords.txt
	cd scripts/dynein/draw && python2 create_MD_array.py


scripts/dynein/draw/tail.py: scripts/dynein/draw/tailDomain.py
	cd scripts/dynein/draw && python2 tailDomain.py




clean:
	rm -f src/*.o
	rm -f src/version-info.h
	rm -f generate_stepping_data
