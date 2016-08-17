import everything, numpy

import matplotlib.pyplot as plt

e = everything.read('5e11_equal_legs')

everything.hist(e.bba, 'bba', e.eqbba)
everything.hist(e.bma, 'bma', e.eqbma)
everything.hist(e.ta, 'ta', e.eqta)
everything.hist(e.uma, 'uma', e.equma)

plt.legend(loc='best')

plt.savefig('plots/5e11_equal_legs_angle_histogram.pdf')
