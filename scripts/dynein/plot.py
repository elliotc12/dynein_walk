"""Functions to make plots for dynein simulation data

NOTE: that there is some weirdness with matplotlib's
pcolor function that chops off some data.
"""



import numpy as np
import matplotlib.pyplot as plt


def getBinIndex(p, bins):
    """Get the index of the bin corresponding to value p
    """

    # deal with edge cases
    if p<= bins[0]:
        return 0
    elif p >= bins[len(bins)-1]:
        return len(bins)-1
    # bin assignment is to the right e.g. if bins are [5, 10] then 7 goes with 10 and 2 goes with 5
    for i in range(0, len(bins)-1):
        if p > bins[i] and p <= bins[i+1]:
            return i
    print("Couldn't find bin for value: {0}".format(p))
    print("Bins: {}".format(bins))
    assert(False) # get cranky if we can't find a bin

def getBinIndex2(p, bins):
    """Return index of bin closest to value p"""

    return np.abs(bins-value).argmin()


def getCounts(X,Y, x_startFromZero=False, y_startFromZero=True, numBins=50):
    """Bin the x and y data"""

    X = np.asarray(X)
    Y = np.asarray(Y)

    print(X.max(), X.min(), Y.max(), Y.min())
    # set up the bins
    if x_startFromZero == True:
        xbins = np.linspace(0, X.max(), numBins+1)
    else:
        xbins = np.linspace(X.min(), X.max(), numBins+1)

    if y_startFromZero == True:
        ybins = np.linspace(0, Y.max(), numBins+1)
    else:
        ybins = np.linspace(Y.min(), Y.max(), numBins+1)

    if len(X) != len(Y):
        print("WARNING: Inputs not of same length--- {0} and {1}".format(len(X), len(Y)))

    counts = np.zeros((len(xbins), len(ybins)))
    for k in range(0, len(X)):
        i = getBinIndex(X[k], xbins)
        j = getBinIndex(Y[k], ybins)
        counts[j,i]+= 1 # rows are y values and columns are x values
    return (xbins, ybins, counts)


def hist_2D(x, y,
            graph_label,
            x_label, y_label,
            colorMap = 'plasma',
            xIsTimeValue=False, yIsTimeValue=False,
            filename=None,
            drawCorrelation=False,
            drawYildizFit=False,
            numBins=50):
    """Create a 2D histogram of x and y data"""

    (xbins, ybins, counts) = getCounts(x, y, xIsTimeValue, yIsTimeValue, numBins)

    fig = plt.figure()
    ax = plt.gca()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(xbins[0], xbins[-1])
    plt.ylim(ybins[0], ybins[-1])
    plt.title(graph_label)

    plt.pcolor(xbins, ybins, counts, cmap=colorMap)
    cb = plt.colorbar()
    cb.set_label('counts')

    if drawCorrelation == True:
        # for explanation of ordinary least squares fitting, check ./papers/notes

        # first create Vandermonde matrix using x values
        A = np.transpose(np.vstack([np.ones(len(x)), x]))
        # solve system A^T*A*b = A^T*y where b is a vector of the parameters
        b0, b1 = np.linalg.lstsq(A, y)[0]
        eq = "Model: y = {:.2} + {:.2}x".format(b0, b1)
        if b1<0:
            eq = "Model: y={:.2} - {:.2}x".format(b0, -b1) # pretty up the case for negative slope
        plt.plot(xbins, b0+b1*xbins, label=eq, linestyle=":", color="white")
        plt.legend()

    # if filename == None:
    #     filename = 'plots/'+graph_label.replace(' ', '_')+".pdf"
    # plt.savefig(filename)

    return fig, ax


if __name__ == "__main__":
    pass 
