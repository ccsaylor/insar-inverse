from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as mtri

"""
Generates plots oof the error of a model as a function of each parameter
"""

filename = 'sensanalysis.txt'
data = np.genfromtxt(filename, dtype = 'float', delimiter='\t')

params = np.asarray(data[:,0]).reshape(6, 100)

error = np.asarray(data[:,1]).reshape(6, 100)

paramnames = ['x', 'y', 'z', 'strike angle', 'dip angle', 'seismic moment']

paramfig, paramax = plt.subplots(2,3, facecolor='w', edgecolor='k', sharex=False, sharey=True)

ix = 0
iy = 0
i = 0

for param in paramnames:
    paramax[ix,iy].plot(params[i], error[i])
    paramax[ix,iy].set_title(param)
    if (iy < 2):
        iy += 1
    else:
        iy = 0
        ix += 1
    i += 1

plt.show()