from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as mtri

filename = 'disptest.txt'
data = np.genfromtxt(filename, dtype = 'float', delimiter='\t')
#points = np.genfromtxt('sourcepoints.txt', dtype = 'float', delimiter='\t')

X = np.asarray(data[:,0])
Y = np.asarray(data[:,1])
Z = np.asarray(data[:,2])

Xf = X[~np.isnan(Z)]
Yf = Y[~np.isnan(Z)]
Zf = Z[~np.isnan(Z)]

minimum = np.nanmin(Z)
maximum = np.nanmax(Z)
Zeros = np.zeros(len(Zf))

fig, ax = plt.subplots(1,1, facecolor='w', edgecolor='k', sharex=True, sharey=True, subplot_kw=dict(projection='3d'))

tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(Xf, Yf, Zeros)

test = ax.plot_trisurf(Xf, Yf, Zf, triangles=tri.triangles, cmap = 'magma', vmin = minimum, vmax = maximum)

ax.set_xlabel('x')
ax.set_ylabel('y')

plt.colorbar(test)

plt.show()