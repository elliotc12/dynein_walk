import everything, numpy

import matplotlib.pyplot as plt

e = everything.read('5e11_equal_legs')

everything.hist(e.bba_PE, 'bba PE')
everything.hist(e.bma_PE, 'bma PE')
everything.hist(e.ta_PE, 'ta PE')
everything.hist(e.uma_PE, 'uma PE')

energies = numpy.linspace(0, 4)
plt.plot(energies, numpy.exp(-energies)/energies**.5/numpy.pi**.5, 'k:', label='Boltzmann')
plt.xlim(xmax=4)

plt.legend(loc='best')

plt.savefig('plots/5e11_equal_legs_PE_histogram.pdf')
