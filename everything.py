
import numpy
import matplotlib.pyplot as plt

class Everything:
    pass

colors = ['r','g','b', 'm', 'c']

def hist(data, name, eqval=None):
    y, bin_edges = numpy.histogram(data, bins=200)
    bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])

    norm = sum(y)*(bin_centers[1] - bin_centers[0])
    plt.errorbar(
        bin_centers,
        y/norm,
        yerr = y**0.5/norm,
        marker = '.',
        color=colors[-1],
        label=name,
    )
    if eqval is not None:
        plt.axvline(eqval, color=colors.pop(), linestyle='dashed')
    else:
        colors.pop()

def read(base):
    n = 0
    def nextn():
        nonlocal n
        n += 1
        return n-1

    elements = Everything()

    data = numpy.loadtxt('data/everything_%s.txt' % base)
    with open('data/everything_%s.txt' % base) as f:
        for line in f:
            if ' T:' in line:
                elements.T = float(line.split(': ')[1])
            if ' kT:' in line:
                elements.kT = float(line.split(': ')[1])
            if 'onebound.bba' in line:
                elements.eqbba = float(line.split(': ')[1])
            if 'onebound.bma' in line:
                elements.eqbma = float(line.split(': ')[1])
            if 'onebound.ta' in line:
                elements.eqta = float(line.split(': ')[1])
            if 'onebound.uma' in line:
                elements.equma = float(line.split(': ')[1])

    elements.time = data[:,nextn()]

    elements.bba_PE = data[:,nextn()]
    elements.bma_PE = data[:,nextn()]
    elements.ta_PE  = data[:,nextn()]
    elements.uma_PE = data[:,nextn()]
    elements.bba = data[:,nextn()]
    elements.bma = data[:,nextn()]
    elements.ta  = data[:,nextn()]
    elements.uma = data[:,nextn()]

    elements.bbx = data[:,nextn()]
    elements.bmx = data[:,nextn()]
    elements.tx  = data[:,nextn()]
    elements.umx = data[:,nextn()]
    elements.ubx = data[:,nextn()]
    elements.bby = data[:,nextn()]
    elements.bmy = data[:,nextn()]
    elements.ty  = data[:,nextn()]
    elements.umy = data[:,nextn()]
    elements.uby = data[:,nextn()]

    elements.f_bbx = data[:,nextn()]
    elements.f_bmx = data[:,nextn()]
    elements.f_tx  = data[:,nextn()]
    elements.f_umx = data[:,nextn()]
    elements.f_ubx = data[:,nextn()]
    elements.f_bby = data[:,nextn()]
    elements.f_bmy = data[:,nextn()]
    elements.f_ty  = data[:,nextn()]
    elements.f_umy = data[:,nextn()]
    elements.f_uby = data[:,nextn()]
    return elements
