from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as mtri

"""
Generates 2D heatmaps of the error of a model that results from varying different pairs of parameters
"""

steps = 50

filename = 'twodsensanalysis.txt'
error = np.genfromtxt(filename, dtype = 'float', delimiter='\t')
errors = error.reshape(36, steps**2)

ranges = [4.0, 4.0, 4.0, np.pi/2.0, np.pi/8.0, 5.0e9]
truevals = [48.5774, 43.6555, -10.4041, 5.25892, 0.523487, 8.29375e+09]
allvals = []
plots = []

for j in range(6):
    truevalue = truevals[j]
    rang = ranges[j]
    stepsize = ranges[j]*2/steps;
    vals = []
    for i in range(steps):
        vals.append(truevalue - rang + stepsize*i);
    allvals.append(vals);

paramnames = ['x', 'y', 'z', 'strike angle', 'dip angle', 'seismic moment']

paramfig, paramax = plt.subplots(6,6, facecolor='w', edgecolor='k', sharex=False, sharey=False)

ix = 0
iy = 0
i = 0

for moderror in errors:
    # if ix != iy:
    temp = paramax[ix,iy].contourf(allvals[ix],allvals[iy], moderror.reshape(steps, steps), 10, cmap = 'magma')
    plots.append(temp)
    paramax[ix,iy].scatter(truevals[ix],truevals[iy], s = 1)
    #paramax[ix,iy].set_title(paramnames[ix] + " & " + paramnames[iy])
    paramax[ix,iy].set_xlabel(paramnames[ix])
    paramax[ix,iy].set_ylabel(paramnames[iy])
    if (iy < 5):
        iy += 1
    else:
        iy = 0
        ix += 1
    i += 1

fig, ax = plt.subplots(3,5, facecolor='w', edgecolor='k', sharex=False, sharey=False)

ax[0,0].contourf(allvals[0],allvals[1], errors[1].reshape(steps, steps), 10, cmap = 'magma')
ax[0,0].scatter(truevals[0],truevals[1], s = 1)
ax[0,0].set_xlabel(paramnames[0])
ax[0,0].set_ylabel(paramnames[1])

ax[0,1].contourf(allvals[0],allvals[2], errors[2].reshape(steps, steps), 10, cmap = 'magma')
ax[0,1].scatter(truevals[0],truevals[2], s = 1)
ax[0,1].set_xlabel(paramnames[0])
ax[0,1].set_ylabel(paramnames[2])

ax[0,2].contourf(allvals[0],allvals[3], errors[3].reshape(steps, steps), 10, cmap = 'magma')
ax[0,2].scatter(truevals[0],truevals[3], s = 1)
ax[0,2].set_xlabel(paramnames[0])
ax[0,2].set_ylabel(paramnames[3])

ax[0,3].contourf(allvals[0],allvals[4], errors[4].reshape(steps, steps), 10, cmap = 'magma')
ax[0,3].scatter(truevals[0],truevals[4], s = 1)
ax[0,3].set_xlabel(paramnames[0])
ax[0,3].set_ylabel(paramnames[4])

ax[0,4].contourf(allvals[0],allvals[5], errors[5].reshape(steps, steps), 10, cmap = 'magma')
ax[0,4].scatter(truevals[0],truevals[5], s = 1)
ax[0,4].set_xlabel(paramnames[0])
ax[0,4].set_ylabel(paramnames[5])

ax[1,0].contourf(allvals[1],allvals[2], errors[8].reshape(steps, steps), 10, cmap = 'magma')
ax[1,0].scatter(truevals[1],truevals[2], s = 1)
ax[1,0].set_xlabel(paramnames[1])
ax[1,0].set_ylabel(paramnames[2])

ax[1,1].contourf(allvals[1],allvals[3], errors[9].reshape(steps, steps), 10, cmap = 'magma')
ax[1,1].scatter(truevals[1],truevals[3], s = 1)
ax[1,1].set_xlabel(paramnames[1])
ax[1,1].set_ylabel(paramnames[3])

ax[1,2].contourf(allvals[1],allvals[4], errors[10].reshape(steps, steps), 10, cmap = 'magma')
ax[1,2].scatter(truevals[1],truevals[4], s = 1)
ax[1,2].set_xlabel(paramnames[1])
ax[1,2].set_ylabel(paramnames[4])

ax[1,3].contourf(allvals[1],allvals[5], errors[11].reshape(steps, steps), 10, cmap = 'magma')
ax[1,3].scatter(truevals[1],truevals[5], s = 1)
ax[1,3].set_xlabel(paramnames[1])
ax[1,3].set_ylabel(paramnames[5])

ax[1,4].contourf(allvals[2],allvals[3], errors[15].reshape(steps, steps), 10, cmap = 'magma')
ax[1,4].scatter(truevals[2],truevals[3], s = 1)
ax[1,4].set_xlabel(paramnames[2])
ax[1,4].set_ylabel(paramnames[3])

ax[2,0].contourf(allvals[2],allvals[4], errors[16].reshape(steps, steps), 10, cmap = 'magma')
ax[2,0].scatter(truevals[2],truevals[4], s = 1)
ax[2,0].set_xlabel(paramnames[2])
ax[2,0].set_ylabel(paramnames[4])

ax[2,1].contourf(allvals[2],allvals[5], errors[17].reshape(steps, steps), 10, cmap = 'magma')
ax[2,1].scatter(truevals[2],truevals[5], s = 1)
ax[2,1].set_xlabel(paramnames[2])
ax[2,1].set_ylabel(paramnames[5])

ax[2,2].contourf(allvals[3],allvals[4], errors[22].reshape(steps, steps), 10, cmap = 'magma')
ax[2,2].scatter(truevals[3],truevals[4], s = 1)
ax[2,2].set_xlabel(paramnames[3])
ax[2,2].set_ylabel(paramnames[4])

ax[2,3].contourf(allvals[3],allvals[5], errors[23].reshape(steps, steps), 10, cmap = 'magma')
ax[2,3].scatter(truevals[3],truevals[5], s = 1)
ax[2,3].set_xlabel(paramnames[3])
ax[2,3].set_ylabel(paramnames[5])

ax[2,4].contourf(allvals[4],allvals[5], errors[29].reshape(steps, steps), 10, cmap = 'magma')
ax[2,4].scatter(truevals[4],truevals[5], s = 1)
ax[2,4].set_xlabel(paramnames[4])
ax[2,4].set_ylabel(paramnames[5])

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)
plt.show()