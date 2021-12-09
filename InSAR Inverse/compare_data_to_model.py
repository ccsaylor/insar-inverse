from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as mtri

"""
Plots 3D displacement maps of the data used for and the model generated by the
InSAR Inverse project
"""

#files = ['mainshock_displacement_binned_50x50_converted.txt', 'mainshock_displacement_binned_50x50_converted_model_data.txt','mainshock_displacement_binned_50x50_converted_residuals.txt']
files = ['mainshock_displacement_binned_converted.txt', 'mainshock_displacement_binned_converted_model_data.txt','mainshock_displacement_binned_converted_residuals.txt']
#files = ['total_displacement_binned_50x50_converted.txt', 'total_displacement_binned_50x50_converted_model_data.txt','total_displacement_binned_50x50_converted_residuals.txt']
#files = ['testing.txt', 'testing_model_data.txt','testing_residuals.txt']
names = ['Data', 'Model','Residuals']
i = 0

# source = np.genfromtxt('converted_data.txt', dtype = 'float', delimiter='\t')
# data = np.asarray(source[:,2])
# minimum = np.amin(data)
# maximum = np.amax(data)

# for filename in files:

    # fig, ax = plt.subplots(1,1, facecolor='w', edgecolor='k', sharex=True, sharey=True, subplot_kw=dict(projection='3d'))

    # sources = np.genfromtxt(filename, dtype = 'float', delimiter='\t')

    # X = np.asarray(sources[:,0])
    # Y = np.asarray(sources[:,1])
    # Z = np.asarray(sources[:,2])
    
    # Xf = X[~np.isnan(Z)]
    # Yf = Y[~np.isnan(Z)]
    # Zf = Z[~np.isnan(Z)]
    # Xcut = [x for x,y,z in zip(Xf, Yf, Zf) if x > 40 and x < 240 and y > 70 and y < 200]
    # Ycut = [y for x,y,z in zip(Xf, Yf, Zf) if x > 40 and x < 240 and y > 70 and y < 200]
    # Zcut = [z for x,y,z in zip(Xf, Yf, Zf) if x > 40 and x < 240 and y > 70 and y < 200]

    # # Zeros = np.zeros(len(Zf))
    # # tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(Xf, Yf, Zeros)
    # # test = ax.plot_trisurf(Xf, Yf, Zf, triangles=tri.triangles, cmap = 'magma')
    
    # Zeros = np.zeros(len(Zcut))
    # tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(Xcut, Ycut, Zeros)
    # test = ax.plot_trisurf(Xcut, Ycut, Zcut, triangles=tri.triangles, cmap = 'magma')

    # ax.set_title(names[i])
    # ax.set_xlabel('x (km)')
    # ax.set_ylabel('y (km)')


    # ax.set_zlim([-0.0006, 0.0006])
    # ax.set_zticklabels([])
    # ax.view_init(azim=-90, elev=90)
    # cbar = plt.colorbar(test)
    # cbar.set_label('Displacement (km)')

    # i += 1

# resfile = np.genfromtxt('mainshock_displacement_binned_converted_residuals.txt', dtype = 'float', delimiter='\t')

# resid = np.asarray(resfile[:,2])
# fig, ax = plt.subplots(1,1, facecolor='w', edgecolor='k', sharex=True, sharey=True)
# #ax.hist(Zf, bins=30)
# ax.hist(resid, bins=30)
# ax.set_xlabel('Residual (km)')
# ax.set_ylabel('Counts')
# plt.show()

class Interferogram:

    def __init__(self, filename):
        self.filename = filename
        self.data = np.genfromtxt(filename, dtype = 'float', delimiter='\t')
        
        X = np.asarray(self.data[:,0])
        Y = np.asarray(self.data[:,1])
        Z = np.asarray(self.data[:,3])*1000
    
        self.X = X[~np.isnan(Z)]
        self.Y = Y[~np.isnan(Z)]
        self.Z = Z[~np.isnan(Z)]
        
    def Subset(self, xmin, xmax, ymin, ymax):
        X = self.X
        Y = self.Y
        Z = self.Z
        self.X = [x for x,y,z in zip(X, Y, Z) if x > xmin and x < xmax and y > ymin and y < ymax]
        self.Y = [y for x,y,z in zip(X, Y, Z) if x > xmin and x < xmax and y > ymin and y < ymax]
        self.Z = [z for x,y,z in zip(X, Y, Z) if x > xmin and x < xmax and y > ymin and y < ymax]
        
class Plotter:
    
    def __init__(self, numx, numy):
        self.numx = numx
        self.numy = numy
        self.fig, self.ax = plt.subplots(numx, numy, facecolor='w', edgecolor='k', subplot_kw=dict(projection='3d'))
        self.fig.set_size_inches(20,12)
        
    def Trisurf(self, plotid, interf, name):
        Zeros = np.zeros(len(interf.Z))
        tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(interf.X, interf.Y, Zeros)
        self.ax[plotid//self.numx, plotid % self.numx].plot_trisurf(interf.X, interf.Y, interf.Z, triangles=tri.triangles, cmap = 'magma')
        self.ax[plotid//self.numx, plotid % self.numx].set_xlabel('x (km)', labelpad=15)
        self.ax[plotid//self.numx, plotid % self.numx].tick_params(pad=5)
        self.ax[plotid//self.numx, plotid % self.numx].set_ylabel('y (km)', labelpad=15)
        self.ax[plotid//self.numx, plotid % self.numx].set_zlabel('z (m)', labelpad=15)
        self.ax[plotid//self.numx, plotid % self.numx].set_title(name)
        self.ax[plotid//self.numx, plotid % self.numx].view_init(azim=45, elev=20)
        #self.ax[plotid//self.numx, plotid % self.numx].set_zlim([-0.5, 0.8])
        for label in self.ax[plotid//self.numx, plotid % self.numx].xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        for label in self.ax[plotid//self.numx, plotid % self.numx].yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        
    def TrisurfTop(self, plotid, interf, name):
        Zeros = np.zeros(len(interf.Z))
        tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(interf.X, interf.Y, Zeros)
        surface = self.ax[plotid//self.numx, plotid % self.numx].plot_trisurf(interf.X, interf.Y, interf.Z, triangles=tri.triangles, cmap = 'magma')
        self.ax[plotid//self.numx, plotid % self.numx].set_xlabel('x (km)', labelpad=32)
        self.ax[plotid//self.numx, plotid % self.numx].tick_params(pad=8)
        self.ax[plotid//self.numx, plotid % self.numx].set_ylabel('y (km)', labelpad=18)
        self.ax[plotid//self.numx, plotid % self.numx].set_title(name)
        self.ax[plotid//self.numx, plotid % self.numx].view_init(azim=-90, elev=89)
        self.ax[plotid//self.numx, plotid % self.numx].set_zticklabels([])
        #self.ax[plotid//self.numx, plotid % self.numx].set_zlim([-0.5, 0.8])
        
        for label in self.ax[plotid//self.numx, plotid % self.numx].xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        for label in self.ax[plotid//self.numx, plotid % self.numx].yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        
        if (plotid == self.numx*self.numy - 1):
            self.fig.colorbar(surface, ax=self.ax.ravel().tolist())
            
class PlotterOneRow:
    
    def __init__(self, num):
        self.num = num
        self.fig, self.ax = plt.subplots(1, num, facecolor='w', edgecolor='k', subplot_kw=dict(projection='3d'))
        self.fig.set_size_inches(20,12)
        
    def Trisurf(self, plotid, interf):
        Zeros = np.zeros(len(interf.Z))
        tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(interf.X, interf.Y, Zeros)
        self.ax[plotid].plot_trisurf(interf.X, interf.Y, interf.Z, triangles=tri.triangles, cmap = 'magma')
        self.ax[plotid].set_xlabel('x (km)', labelpad=15)
        self.ax[plotid].tick_params(pad=10)
        self.ax[plotid].set_ylabel('y (km)', labelpad=15)
        self.ax[plotid].set_zlabel('z (m)', labelpad=15)
        self.ax[plotid].view_init(azim=45, elev=20)
        #self.ax[plotid].set_zlim([-0.5, 0.8])
        for label in self.ax[plotid].xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        for label in self.ax[plotid].yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        
    def TrisurfTop(self, plotid, interf):
        Zeros = np.zeros(len(interf.Z))
        tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(interf.X, interf.Y, Zeros)
        surface = self.ax[plotid].plot_trisurf(interf.X, interf.Y, interf.Z, triangles=tri.triangles, cmap = 'magma')
        self.ax[plotid].set_xlabel('x (km)', labelpad=20)
        self.ax[plotid].tick_params(pad=12)
        self.ax[plotid].set_ylabel('y (km)', labelpad=35)
        self.ax[plotid].view_init(azim=-90, elev=89)
        self.ax[plotid].set_zticklabels([])
        #self.ax[plotid].set_zlim([-0.5, 0.8])
        for label in self.ax[plotid].xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        for label in self.ax[plotid].yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        
        if (plotid == self.num - 1):
            self.fig.colorbar(surface, ax=self.ax.ravel().tolist())

if __name__ == "__main__":

    plt.rcParams['font.size'] = '18'

    intdata = Interferogram('iran_binned_fitted_data.txt')
    #intdata.Subset(553, 620, 3837, 3880)
    
    intmodel = Interferogram('iran_binned_model_data.txt')
    #intmodel.Subset(553, 620, 3837, 3880)
    
    plotter = Plotter(2, 2)
    plotter.Trisurf(0, intdata, "Data")
    plotter.TrisurfTop(1, intdata, "Data")
    plotter.Trisurf(2, intmodel, "Model")
    plotter.TrisurfTop(3, intmodel, "Model")
    
    # intresid = Interferogram('mainshock_displacement_binned_converted_residuals.txt')
    # intresid.Subset(40, 240, 70, 200)
    
    # resplotter = PlotterOneRow(2)
    # resplotter.Trisurf(0, intresid)
    # resplotter.TrisurfTop(1, intresid)

    plt.show()
