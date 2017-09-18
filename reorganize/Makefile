CPPFLAGS = -std=c++98 -g -Werror -O2 -Wall
LIBRARIES = -lm

.PHONY: clean

.PRECIOUS: data/stepping_data_%.txt data/stepping_config_%.txt data/stepping_movie_data_%.txt data/bothbound_data_%.bin data/onebound_data_%.bin data/ob_config_%.txt data/bb_config_%.txt # prevent nonexistant data files from being deleted after creation+use

all: generate_stepping_data

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

generate_stepping_data: src/generate_stepping_data.cpp src/dynein_simulate.o src/dynein_struct_onebound.o src/dynein_struct_bothbound.o src/utilities.o src/version-info.h
	mkdir -p runlogs
	mkdir -p data
	$(CXX) -c -o src/generate_stepping_data.o src/generate_stepping_data.cpp $(CPPFLAGS)
	$(CXX) src/generate_stepping_data.o src/dynein_simulate.o src/dynein_struct_onebound.o src/dynein_struct_bothbound.o src/utilities.o -o generate_stepping_data

clean:
	rm -f *.o src/*.o
	rm -f generate_stepping_data
