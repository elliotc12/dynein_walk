# dynein_walk

[![Build Status](https://semaphoreci.com/api/v1/droundy/dynein_walk/branches/master/shields_badge.svg)](https://semaphoreci.com/droundy/dynein_walk)
[![build status](https://gitlab.com/daveroundy/dynein_walk/badges/master/build.svg)](https://gitlab.com/daveroundy/dynein_walk/commits/master)

A simulation of the stepping pattern of the molecular motor Dynein using Monte Carlo methods and Brownian mechanics.

## Building

To build this project you will need
[`fac`](http://physics.oregonstate.edu/~roundyd/fac/), inkscape,
texlive, (specifically, bibtex and pdflatex), python3 with matplotlib,
and a C++ compiler such as g++.  On a Debian-derived system you can
install all of these but `fac` with:

    apt-get install inkscape texlive texlive-xetex texlive-latex-extra build-essential python3-matplotlib
