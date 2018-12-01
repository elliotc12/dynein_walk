import glob

for svg in glob.glob('*.svg'):
    basename = svg[:-4]
    print('| inkscape -D --export-pdf {basename}.pdf {basename}.svg'.format(basename=basename))
