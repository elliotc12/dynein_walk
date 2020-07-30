import glob

for svg in glob.glob('*.svg'):
    basename = svg[:-4]
    print('| cairosvg -o {basename}.pdf {basename}.svg'.format(basename=basename))
