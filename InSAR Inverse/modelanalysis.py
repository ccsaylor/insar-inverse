from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as mtri

"""
Generates histograms of each parameter resulting from different runs of a genetic algorithm
fit
"""

filenamenoext = 'interferogram'
data = np.genfromtxt(filenamenoext + '_fitted_data.txt', dtype = 'float', delimiter='\t')
datasources = np.genfromtxt(filenamenoext + '_sources.txt', dtype = 'float', delimiter='\t')
modeldata = np.genfromtxt(filenamenoext + '_modeldata.txt', dtype = 'float', delimiter='\t')
modelsources = np.genfromtxt(filenamenoext + '_modelsources.txt', dtype = 'float', delimiter='\t')
errors = np.genfromtxt(filenamenoext + '_errors.txt', dtype = 'float', delimiter='\t')

X = np.asarray(data[:,0])
Y = np.asarray(data[:,1])
Z = np.asarray(data[:,2])

truevals = [48.5774, 43.6555, -10.4041, 5.25892, 0.523487, 8.29375e+09]

cutmodeldata = np.asarray([mdata for mdata,msources,merrors in zip(modeldata, modelsources, errors) if merrors < 0.0001])
cutmodelsources = np.asarray([msources for mdata,msources,merrors in zip(modeldata, modelsources, errors) if merrors < 0.0001])
cuterrors = np.asarray([merrors for mdata,msources,merrors in zip(modeldata, modelsources, errors) if merrors < 0.0001])

print(len(cutmodeldata))
modeldata = cutmodeldata
modelsources = cutmodelsources
errors = cuterrors

Zeros = np.zeros(len(Z))

# fig, ax = plt.subplots(1,2, facecolor='w', edgecolor='k', sharex=True, sharey=True, subplot_kw=dict(projection='3d'))

# tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(X, Y, Zeros)
# test = ax[0].plot_trisurf(X, Y, Z, triangles=tri.triangles, cmap = 'magma')
# ax[1].plot_trisurf(X, Y, modeldata[1], triangles=tri.triangles, cmap = 'magma')

# plt.colorbar(test)

# for axis in ax:
    # axis.set_xlabel('x')
    # axis.set_ylabel('y')
    # axis.set_zticklabels([])
    # axis.view_init(azim=-90, elev=90)

# plt.show()


params = ['x', 'y', 'z', 'strike angle', 'dip angle', 'seismic moment']

num_models = len(errors)

indicesup = []
indicesdown = []

for k in range(len(modelsources[:,3])):
    if modelsources[k,3] < 0:
        indicesup.append(k)
    elif modelsources[k,3] > 2*np.pi:
        indicesdown.append(k)
for index in indicesup:
    modelsources[index][3] += 2*np.pi

for index in indicesdown:
    modelsources[index][3] -= 2*np.pi

paramfig, paramax = plt.subplots(2,3, facecolor='w', edgecolor='k', sharex=False, sharey=False)

ix = 0
iy = 0
i = 0
for param in params:
    paramax[ix,iy].scatter(modelsources[:,i], np.zeros(num_models))
    paramax[ix,iy].set_title(param)
    if (iy < 2):
        iy += 1
    else:
        iy = 0
        ix += 1
    i += 1

ix = 0
iy = 0
i = 0
histfig, histax = plt.subplots(2,3, facecolor='w', edgecolor='k', sharex=False, sharey=False)

for param in params:
    histax[ix,iy].hist(modelsources[:,i])
    histax[ix,iy].set_title(param)
    if i != 5:
        histax[ix,iy].text(0.8, 0.9, 'σ = %.3f' % np.std(modelsources[:,i]), horizontalalignment='center', verticalalignment='center', transform=histax[ix,iy].transAxes)
        histax[ix,iy].text(0.8, 0.8, 'Mean = %.3f' % np.mean(modelsources[:,i]), horizontalalignment='center', verticalalignment='center', transform=histax[ix,iy].transAxes)
        histax[ix,iy].text(0.8, 0.7, 'True = %.3f' % truevals[i], horizontalalignment='center', verticalalignment='center', transform=histax[ix,iy].transAxes)
    if (iy < 2):
        iy += 1
    else:
        iy = 0
        ix += 1
    i += 1

histax[1,2].text(0.8, 0.9, 'σ = ' + "{:.2e}".format(np.std(modelsources[:,5])), horizontalalignment='center', verticalalignment='center', transform=histax[1,2].transAxes)
histax[1,2].text(0.8, 0.8, 'Mean = ' + "{:.2e}".format(np.mean(modelsources[:,5])), horizontalalignment='center', verticalalignment='center', transform=histax[1,2].transAxes)
histax[1,2].text(0.8, 0.7, 'True = ' + "{:.2e}".format(truevals[5]), horizontalalignment='center', verticalalignment='center', transform=histax[1,2].transAxes)
plt.show()