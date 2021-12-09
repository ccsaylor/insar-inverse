from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import matplotlib.tri as mtri
from sklearn import linear_model
from lmfit import Minimizer, Parameters, fit_report
import time
from scipy.interpolate import RegularGridInterpolator
    
class DataStore:
    def __init__(self):
        self.test = []
    
    def ImportData(self, filename, xmin = -np.inf, xmax = np.inf, ymin = -np.inf, ymax = np.inf):
        self.data = np.genfromtxt(filename, dtype = 'float', delimiter='\t')
        self.y = np.concatenate((self.data[:,2], self.data[:,3]), axis=None)
        self.yew = self.data[:,2]
        self.yud = self.data[:,3]
        self.numdata = len(self.data)
        print(self.numdata, "data points.")
        
    def ImportPoints(self, filename):
        points = np.genfromtxt(filename, dtype = 'float', delimiter='\t')
        
        self.numpoints = len(points)
        self.px = np.asarray(points[:,0])
        self.py = np.asarray(points[:,1])
        self.pz = np.asarray(points[:,2])
        self.pstrike = np.asarray(points[:,3])
        self.pdip = np.asarray(points[:,4])
        self.prake = np.asarray(points[:,5])
        self.pm = np.asarray(points[:,6])
        
        print(self.numpoints, "sources.")
        
    def CalculateGreens(self, filename):
        greens = np.genfromtxt(filename, dtype='float', delimiter='\t')
        x = greens[:,0]
        y = greens[:,1]
        ewstrike = greens[:,2]
        ewdip = greens[:,3]
        udstrike = greens[:,4]
        uddip = greens[:,5]
        
        ux = np.unique(x)
        uy = np.unique(y)
        
        ewstrint = RegularGridInterpolator((ux, uy), ewstrike.reshape(len(uy), len(ux)), bounds_error=False, fill_value=0)
        ewdipint = RegularGridInterpolator((ux, uy), ewdip.reshape(len(uy), len(ux)), bounds_error=False, fill_value=0)
        udstrint = RegularGridInterpolator((ux, uy), udstrike.reshape(len(uy), len(ux)), bounds_error=False, fill_value=0)
        uddipint = RegularGridInterpolator((ux, uy), uddip.reshape(len(uy), len(ux)), bounds_error=False, fill_value=0)
        
        self.ewstrgreens = []
        self.ewdipgreens = []
        self.udstrgreens = []
        self.uddipgreens = []
        for dx, dy in zip(self.data[:,0], self.data[:,1]):
        
            newx = []
            newy = []
            for px, py, pstr in zip(self.px, self.py, self.pstrike):
                xtemp = (dx - px)*np.cos(pstr) - (dy - py)*np.sin(pstr)
                ytemp = (dx - px)*np.sin(pstr) + (dy - py)*np.cos(pstr)
                newx.append(xtemp)
                newy.append(ytemp)
        
            newpoints = np.stack((newx, newy), axis=-1)
            ewstrtemp = ewstrint(newpoints)
            ewdiptemp = ewdipint(newpoints)
            udstrtemp = udstrint(newpoints)
            uddiptemp = uddipint(newpoints)
            
            self.ewstrgreens.append(ewstrtemp)
            self.ewdipgreens.append(ewdiptemp)
            self.udstrgreens.append(udstrtemp)
            self.uddipgreens.append(uddiptemp)
            
        self.ewstrgreens = np.asarray(self.ewstrgreens)
        self.ewdipgreens = np.asarray(self.ewdipgreens)
        self.udstrgreens = np.asarray(self.udstrgreens)
        self.uddipgreens = np.asarray(self.uddipgreens)
        
    def ImportBeta(self, filename):
        points = np.genfromtxt(filename, dtype = 'float', delimiter='\t')
        moments = np.asarray(points[:,6])
        rakes = np.asarray(points[:,5])
        self.Beta = np.concatenate((moments, rakes), axis=0)
        
    def fnc2min(self, params):
        coeffs = []
        for key, value in params.valuesdict().items():
            coeffs.append(value)
            
        moments = coeffs[:self.numpoints]
        rakes = coeffs[self.numpoints:]
            
        ewmodel = []
        udmodel = []

        for i in range(len(self.yew)):
            
            ewstrgreens = self.ewstrgreens[i,:]
            ewdipgreens = self.ewdipgreens[i,:]
            udstrgreens = self.udstrgreens[i,:]
            uddipgreens = self.uddipgreens[i,:]
            
            ewgreens = np.multiply(ewstrgreens, np.cos(rakes)) + np.multiply(ewdipgreens, np.sin(rakes))
            udgreens = np.multiply(udstrgreens, np.cos(rakes)) + np.multiply(uddipgreens, np.sin(rakes))
            
            ewmodel.append(np.sum(np.multiply(moments, ewgreens)))
            udmodel.append(np.sum(np.multiply(moments, udgreens)))

        return ewmodel - self.yew + udmodel - self.yud

    def CalculateBetaLM(self, max=np.inf):
    
        indices = np.arange(0, len(self.yew))
        
        params = Parameters()
        for j in range(self.numpoints):
            name = "p%i" % j
            val = np.random.uniform(1.0e6, 1.0e9)
            params.add(name, value=val, min=0, max = max)
            
        for j in range(self.ewstrgreens.shape[1]):
            name = "p%i" % j
            params.add(name + "rake", value=self.prake[j], vary=False)
            
        minner = Minimizer(self.fnc2min, params, fcn_args=())
        result = minner.minimize()
        
        self.Beta = np.asarray([value for key, value in result.params.valuesdict().items()])

    def SaveModel(self, outfilename):
        
        with open(outfilename, 'w') as ofile:
            for px, py, pz, pstrike, pdip, prake, pm in zip(self.px, self.py, self.pz, self.pstrike, self.pdip, self.prake, self.Beta):
                ofile.write(str(px) + "\t" + str(py) + "\t" + str(pz) + "\t" + str(pstrike) + "\t" + str(pdip) + "\t" + str(prake) + "\t" + str(pm) + "\n")

class Plotter:
    
    def __init__(self, store, num_plots):
        self.datastore = store
        self.fig, self.ax = plt.subplots(1, num_plots, facecolor='w', edgecolor='k', subplot_kw=dict(projection='3d'))
        self.dataid = -1
        self.vmin = 0
        self.vmax = 0

    def PlotData(self, plotid):
    
        self.dataid = plotid
    
        Zeros = np.zeros(len(self.datastore.data))
        tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(self.datastore.data[:,0], self.datastore.data[:,1], Zeros)
        
        self.ax[plotid].plot_trisurf(np.add(self.datastore.data[:,0], self.datastore.data[:,2]), self.datastore.data[:,1], self.datastore.data[:,3], triangles=tri.triangles, cmap = 'magma')
        self.vmin = self.datastore.data[:,3].min()
        self.vmax = self.datastore.data[:,3].max()
        self.ax[plotid].set_xlabel('x (km)')
        self.ax[plotid].set_ylabel('y (km)')
        self.ax[plotid].set_zlabel('z (m)')
        self.ax[plotid].set_title("Data")
        #self.ax[plotid].view_init(azim=45, elev=20)
        self.ax[plotid].view_init(azim=-90, elev=90)
        
    def PlotModel(self, plotid):
    
        Zeros = np.zeros(len(self.datastore.data))
        tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(self.datastore.data[:,0], self.datastore.data[:,1], Zeros)
        coeffs = self.datastore.Beta
        moments = coeffs[:self.datastore.numpoints]
        rakes = coeffs[self.datastore.numpoints:]
        
        ewmodel = []
        udmodel = []
        
        for i in range(len(self.datastore.yew)):
            
            ewstrgreens = self.datastore.ewstrgreens[i,:]
            ewdipgreens = self.datastore.ewdipgreens[i,:]
            udstrgreens = self.datastore.udstrgreens[i,:]
            uddipgreens = self.datastore.uddipgreens[i,:]
            
            ewgreens = np.multiply(ewstrgreens, np.cos(rakes)) + np.multiply(ewdipgreens, np.sin(rakes))
            udgreens = np.multiply(udstrgreens, np.cos(rakes)) + np.multiply(uddipgreens, np.sin(rakes))
            
            ewmodel.append(np.sum(np.multiply(moments, ewgreens)))
            udmodel.append(np.sum(np.multiply(moments, udgreens)))
    
        self.ax[plotid].plot_trisurf(np.add(self.datastore.data[:,0], ewmodel), self.datastore.data[:,1], udmodel, triangles=tri.triangles, cmap = 'magma', vmin = self.vmin, vmax = self.vmax)
        self.ax[plotid].set_xlabel('x (km)')
        self.ax[plotid].set_ylabel('y (km)')
        self.ax[plotid].set_zlabel('z (m)')
        self.ax[plotid].set_title("Model")
        #self.ax[plotid].view_init(azim=45, elev=20)
        self.ax[plotid].view_init(azim=-90, elev=90)
        
    def PlotModelPoints(self, plotid):
    
        moments = self.datastore.Beta[:self.datastore.numpoints]
    
        px = self.datastore.px[np.where(moments != 0)]
        py = self.datastore.py[np.where(moments != 0)]
        pz = self.datastore.pz[np.where(moments != 0)]
        m = moments[np.where(moments != 0)]
    
        scatter = self.ax[plotid].scatter(px, py, pz, c = m, s = 20)

        self.ax[plotid].set_xlabel('x (km)')
        self.ax[plotid].set_ylabel('y (km)')
        self.ax[plotid].set_zlabel('z (m)')
        #self.ax[plotid].view_init(azim=45, elev=20)
        self.ax[plotid].view_init(azim=-90, elev=90)
        
        if self.dataid != -1:
            self.ax[plotid].set_xlim(self.ax[self.dataid].get_xlim())
            self.ax[plotid].set_ylim(self.ax[self.dataid].get_ylim())
        #plt.colorbar(scatter)
        
if __name__ == "__main__":

    start = time.time()
    
    store = DataStore()
    # store.ImportData("ridgecrest_displacements_100x100_fitted_data.txt")
    # store.ImportPoints("ridgecrest_displacements_100x100_model_points.txt")
    store.ImportData("iran_displacements_100x100_fitted_data.txt")
    store.ImportPoints("iran_displacements_100x100_model_points.txt")
    store.CalculateGreens("disptest.txt")
    
    store.CalculateBetaLM()
    store.SaveModel("2doptimizedmodel.txt")
    
    plotter = Plotter(store, 3)
    plotter.PlotData(0)
    plotter.PlotModel(1)
    plotter.PlotModelPoints(2)
    
    plt.show()
    
    end = time.time()
    print(end - start)
    
