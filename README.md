# dynein_walk

[![Build Status](https://semaphoreci.com/api/v1/droundy/dynein_walk/branches/master/shields_badge.svg)](https://semaphoreci.com/droundy/dynein_walk)
[![build status](https://gitlab.com/daveroundy/dynein_walk/badges/master/build.svg)](https://gitlab.com/daveroundy/dynein_walk/commits/master)

A simulation of the stepping pattern of the molecular motor Dynein using Brownian mechanics.

## Summer 2019 to-do items

- Research
  - Try to get `bb-mc.py` actually working, which probably means finding a suspected bug in `onebound.cpp`.
  - Look through how `onebound.cpp` is working internally.
  - Create figures (rinse and repeat)
  - Create movie data (and movies)?
  - Generate time simulation (i.e. string together probabilities in MC sequence)

- Write latex description of method and figures (eventually thesis)

- Read dynein papers & summarize in the thesis
  - Annotated bibliography
  - Bibtex file
  - Google scholar

- Obtain data from papers
  - Email authors
  - Use gimp or similar? :(

- Read code
  - `onebound.cpp`
  - `dynein_onebound_struct.cpp`
  - *add more later*

- Check linear algebra from Mathematica for the onebound solution.
  - backward may be easier
  - forward may give a nicer result
  - consider using different set of angles (e.g. relative angles)

- Reimplement `onebound` (in rust?)
  - Note that this may be a way to track down a bug in the C++ code.

- Implement toy dynamics/toy MC for learning

## Building

To build this project you will need
[`fac`](http://physics.oregonstate.edu/~roundyd/fac/), inkscape,
texlive, (specifically, bibtex and pdflatex), python3 with matplotlib,
and a C++ compiler such as g++.  On a Debian-derived system you can
install all of these but `fac` with:

    apt-get install inkscape texlive texlive-xetex texlive-latex-extra build-essential python3-matplotlib
