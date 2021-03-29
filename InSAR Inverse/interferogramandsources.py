from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as mtri

fig, ax = plt.subplots(1,1, facecolor='w', edgecolor='k', sharex=True, sharey=True, subplot_kw=dict(projection='3d'))
  
interferogram = np.genfromtxt('interf1test.txt', dtype = 'float', delimiter='\t')
points = np.genfromtxt('interf1test_sourcepoints.txt', dtype = 'float', delimiter='\t')

x = np.asarray(interferogram[:,0])
y = np.asarray(interferogram[:,1])
z = np.asarray(interferogram[:,2])
zeros = np.zeros(len(z))

# a = np.asarray(points[:,0]) #use for multiple point sources
# b = np.asarray(points[:,1])
# c = np.asarray(points[:,2])

a = np.asarray(points[0]) #use for single point source
b = np.asarray(points[1])
c = np.asarray(points[2])

tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(x, y, zeros)
ax.plot_trisurf(np.asarray(x), np.asarray(y), np.asarray(z), triangles=tri.triangles, cmap = 'magma') 
ax.scatter(a, b, c)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()