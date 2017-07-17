from __future__ import division
import numpy as np
import sys

''' 
This program will look into the specified file and create a new python file with arrays for the data columns: times, PE_1, PE_2, PE_3, PE_4, and PE_5

example for running script from terminal: python autocorrelation.py data.txt #
expecting thesis_movie_data.txt
'''

path = sys.argv[1]
dataTable = np.loadtxt(path, delimiter='\t', skiprows=1)
print "Successfully loaded!"
#print np.shape(dataTable) --- 

times = dataTable[:,1] 
PE_1 = dataTable[:,2]
PE_2 = dataTable[:,3]
PE_3 = dataTable[:,4]
PE_4 = dataTable[:,5]
PE_5 = dataTable[:,6]

print repr(PE_1)

with open("data_arrays.py" , "w") as f:
    f.write("""from numpy import array \n""")
    f.write("times = %s \n" % repr(times))
    f.write("PE_1 = %s \n" % repr(PE_1))
    f.write("PE_2 = %s \n" % repr(PE_2))
    f.write("PE_3 = %s \n" % repr(PE_3))
    f.write("PE_4 = %s \n" % repr(PE_4))
    f.write("PE_5 = %s \n" % repr(PE_5))

print "All done!" 
