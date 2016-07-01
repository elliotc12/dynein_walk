#!/usr/bin/python2

import glob

cppfiles = glob.glob('*.cpp') + glob.glob('simulations/*.cpp')
mainfiles = []
for cpp in cppfiles:
    print '| g++ -std=c++11 -g -Wall -Werror -O2 -c -o %s %s' % (cpp[:-3]+'o', cpp)
    print '> %s.o' % cpp[:-4]
    print
    with open(cpp) as f:
        if f.read().find('main(') >= 0:
            mainfiles.append(cpp)

otherfiles = [f for f in cppfiles if f not in mainfiles]
ofiles = [f[:-3]+'o' for f in otherfiles]

for main in mainfiles:
    exename = main[:-4]
    if exename.find('/') > 0:
        exename = exename[exename.find('/')+1:]
    print '| g++ %s %s -o %s' % (main[:-3]+'o', ' '.join(ofiles), exename)
    print '>', exename
    for c in cppfiles:
        print '<', c[:-3]+'o'
    print
