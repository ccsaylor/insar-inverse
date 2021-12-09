from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as mtri
from scipy.interpolate import RegularGridInterpolator

class Interpolator:

    def __init__(self, filename):
    
        data = np.genfromtxt(filename, dtype = 'float', delimiter='\t')

        x = np.asarray(data[:,0])
        y = np.asarray(data[:,1])
        ewstr = np.asarray(data[:,2])
        ewdip = np.asarray(data[:,3])
        udstr = np.asarray(data[:,4])
        uddip = np.asarray(data[:,5])

        ux = np.unique(x)
        uy = np.unique(y)

        self.ewstrint = RegularGridInterpolator((ux, uy), ewstr.reshape(len(uy), len(ux)), bounds_error=False, fill_value=0)
        self.ewdipint = RegularGridInterpolator((ux, uy), ewdip.reshape(len(uy), len(ux)), bounds_error=False, fill_value=0)
        self.udstrint = RegularGridInterpolator((ux, uy), udstr.reshape(len(uy), len(ux)), bounds_error=False, fill_value=0)
        self.uddipint = RegularGridInterpolator((ux, uy), uddip.reshape(len(uy), len(ux)), bounds_error=False, fill_value=0)
        
class Source:
    
    def __init__(self, x, y, strike, rake, moment):
    
        self.x = x
        self.y = y
        self.strike = strike
        self.rake = rake
        self.m = moment

    def CalcDisplacement(self, interp, xs, ys):
    
        newx = []
        newy = []
        for x, y in zip(xs, ys):
            newx.append((x - self.x)*np.cos(self.strike) - (y - self.y)*np.sin(self.strike))
            newy.append((x - self.x)*np.sin(self.strike) + (y - self.y)*np.cos(self.strike))
        newpoints = np.stack((newx, newy), axis=-1)
        
        ewstrinterp = interp.ewstrint(newpoints)
        ewdipinterp = interp.ewdipint(newpoints)
        udstrinterp = interp.udstrint(newpoints)
        uddipinterp = interp.uddipint(newpoints)
        
        ewdisp = ewstrinterp*np.cos(self.rake) + ewdipinterp*np.sin(self.rake)
        uddisp = udstrinterp*np.cos(self.rake) + uddipinterp*np.sin(self.rake)
        
        return ewdisp, uddisp
            
        
if __name__ == "__main__":

    inter = Interpolator("single.txt")

    xs = np.linspace(-50, 100, 50)
    ys = np.linspace(-50, 50, 50)
    
    X, Y = np.meshgrid(xs, ys)
    xs = X.flatten()
    ys = Y.flatten()
    
    src = Source(0, 0, 4.0492, 3.2638, 1)
    ewdisp, uddisp = src.CalcDisplacement(inter, xs, ys)
    
    src2 = Source(50, 0, 2.4609, -0.1571, 0)
    ewdisp2, uddisp2 = src2.CalcDisplacement(inter, xs, ys)
    
    ewdisp = ewdisp + ewdisp2
    uddisp = uddisp + uddisp2

    zeros = np.zeros(len(xs))

    fig, ax = plt.subplots(1,1, facecolor='w', edgecolor='k', sharex=True, sharey=True, subplot_kw=dict(projection='3d'))

    tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(xs, ys, zeros)

    test = ax.plot_trisurf(xs + ewdisp, ys, uddisp, triangles=tri.triangles, cmap = 'magma')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.view_init(azim=-90, elev=90)

    plt.show()


