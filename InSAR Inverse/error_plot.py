from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as mtri

"""
Plots the average error of the population of solutions for the
InSAR Inverse project
"""

filename = 'total_displacement_binned_converted_error.txt' #format is generation number -- min error -- average error -- max error separated by tabs

file = open(filename, 'r')

sources = np.genfromtxt(filename, dtype = 'float', delimiter='\t')

generation = np.asarray(sources[:,0])
min_error = np.asarray(sources[:,1])
average_error = np.asarray(sources[:,2])
max_error = np.asarray(sources[:,3])

fig = plt.figure()

#ax1 = fig.add_subplot(131)
#ax1.plot(generation, min_error)

#ax2 = fig.add_subplot(132)
ax2 = fig.add_subplot(111)
ax2.plot(generation, average_error)

#ax3 = fig.add_subplot(133)
#ax3.plot(generation, max_error)

#ax2.set_yscale('log')

ax2.set_xlabel('Generation')
ax2.set_ylabel(r'$\chi^2$')

#plt.savefig('average_error.svg',dpi=400,transparent=True)
plt.show()
