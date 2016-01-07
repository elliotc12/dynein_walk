#!/usr/bin/python2

import glob

cppfiles = glob.glob('*.cpp')
for cpp in cppfiles:
    print '| g++ -std=c++11 -g -Wall -Werror -O2 -c %s' % cpp
    print '> %s.o' % cpp[:-4]
    print

mainfiles = ['dynein_test.cpp', 'dynein_walk.cpp']
otherfiles = [f for f in cppfiles if f not in mainfiles]
ofiles = [f[:-3]+'o' for f in otherfiles]

for main in mainfiles:
    print '| g++ %s %s -o %s' % (main[:-3]+'o', ' '.join(ofiles), main[7:-4])
    print
