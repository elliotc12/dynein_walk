from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


path = r"thesis_movie_data.txt"

dataTable = np.loadtxt(path, skiprows=2, comments='#', delimiter='\t')

print np.shape(dataTable) 
