"""Functions to make plots for dynein simulation data

NOTE: Talk with david about how we are binning... (from left, right, or middle?)
"""



import numpy as np
import matplotlib.pyplot as plt




#--------------- 2 dim Histogram ---------------------------#


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


def hist_2d(x, y,
            colorMap = 'plasma',
            xIsTimeValue=False, yIsTimeValue=False,
            drawCorrelation=False,
            drawYildizFit=False,
            numBins=50):
    """Create a 2D histogram of x and y data

    Returns: fig, ax"""

    (xbins, ybins, counts) = getCounts(x, y, xIsTimeValue, yIsTimeValue, numBins)

    fig = plt.figure()
    ax = plt.gca()
    plt.xlim(xbins[0], xbins[-1])
    plt.ylim(ybins[0], ybins[-1])
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

    return fig, ax


#-------------- trajectory projection plots -----------------------------#
import os, sys
sys.path.insert(0, '../scripts') # need to do this to make import work
import dynein.draw.cartoon as cartoon
from matplotlib import gridspec
from matplotlib.patches import Rectangle

def trajectory(movieData):
    """Create trajectory projection plots of x and y data

    RETURNS: fig, ax
    """
    # note 1 and 5 correspond to the binding domains e.g. x1, y1, x5, y5, etc...
    num_to_avg = 300
    points = len(movieData.times) // num_to_avg # integer division


    # NOTE as of right now, the ordering of domains in the movie data files depends
    # on the state dynein is in e.g. NEARBOUND or FARBOUND as defined in dynein_struct.h


    # average data at each point (to the right). This is done using slices across a multiple of num_to_avg
    # if we just want to get rid of noise we might consider using a convolution with a filter to smooth
    # he data out.
    avg_x1 = np.array([np.mean(movieData.x1[num_to_avg*i:num_to_avg*(i+1)]) for i in range(num_to_avg)])
    avg_x5 = np.array([np.mean(movieData.x5[num_to_avg*i:num_to_avg*(i+1)]) for i in range(num_to_avg)])
    avg_y1 = np.array([np.mean(movieData.y1[num_to_avg*i:num_to_avg*(i+1)]) for i in range(num_to_avg)])
    avg_y5 = np.array([np.mean(movieData.y5[num_to_avg*i:num_to_avg*(i+1)]) for i in range(num_to_avg)])
    avg_times = np.array([np.mean(movieData.times[num_to_avg*i:num_to_avg*(i+1)]) for i in range(num_to_avg)])

    # find global maxima to set plot limits
    max_x = np.max([np.max(avg_x1), np.max(avg_x5)])
    max_y = np.max([np.max(avg_y1), np.max(avg_y5)])
    min_x = np.min([np.min(avg_x1), np.min(avg_x5)])
    min_y = np.min([np.min(avg_y1), np.min(avg_y5)])

    fig = plt.figure()
    # use gridspec to control how subplots are placed
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)

    # turn off ticklabels as we will use the same times for all three graphs
    plt.setp([ax0.get_xticklabels(), ax1.get_xticklabels()], visible=False)

    # x-projection
    ax0.set_ylabel("x-projection (nm)")
    ax0.set_ylim(min_x-1, max_x+1)
    ax0.plot(avg_times, avg_x1, label="near foot", c='b')
    ax0.plot(avg_times, avg_x5, label="far foot", c='r')
    ax0.legend(loc="upper right")

    # NOTE: Need to look closer at how cartoons are generated. It looks like we are currently just choosing

    # y-projection
    ax1.set_ylabel("y-projection (nm)")
    ax1.set_xlabel("time (s)")
    ax1.set_ylim(min_y-1, max_y+1)
    ax1.plot(avg_times, avg_y1, label="near foot", c='b')
    ax1.plot(avg_times, avg_y5, label="far foot", c='r')
    ax1.legend(loc="upper right")
    gs.tight_layout(fig, h_pad=0)

    ax = plt.gca()
    return fig, ax


if __name__ == "__main__":

    import dynein.data as Data

    testData = Data.MovieData("../data/paper_trajectory_movie_data.txt")
    fig, ax = trajectory(testData)

    plt.figure()
    plt.plot(testData.times, testData.x1)
    plt.show()

