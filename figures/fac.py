#!/usr/bin/python2

import glob

for svg in glob.glob('*.svg'):
    print "| inkscape -D --export-pdf %s %s" % (svg[:-3]+'pdf', svg)
    print "C ~/.config"
    print
