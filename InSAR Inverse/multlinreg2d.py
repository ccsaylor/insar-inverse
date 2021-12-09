from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import matplotlib.tri as mtri
from sklearn import linear_model
from lmfit import Minimizer, Parameters, fit_report
import time
from scipy.ndimage import gaussian_filter
    
class DataStore:
    def __init__(self):
        self.test = []
    
    def ImportData(self, filename, xmin = -np.inf, xmax = np.inf, ymin = -np.inf, ymax = np.inf):
        self.data = np.genfromtxt(filename, dtype = 'float', delimiter='\t')
        self.y = np.concatenate((self.data[:,2], self.data[:,3]), axis=None)
        
    def SmoothData(self):
        self.datasmooth = self.data.copy()
        self.datasmooth[:,2] = gaussian_filter(self.datasmooth[:,2], sigma=1)
        self.datasmooth[:,3] = gaussian_filter(self.datasmooth[:,3], sigma=1)
        
        self.ysmooth = np.concatenate((self.datasmooth[:,2], self.datasmooth[:,3]), axis=None)
        
    def ImportGreens(self, eastwest = "greens_x.txt", updown = "greens_z.txt"):
        ew = np.genfromtxt(eastwest, dtype = 'float', delimiter='\t')
        ud = np.genfromtxt(updown, dtype = 'float', delimiter='\t')

        self.X = np.concatenate((ew, ud), axis=0)
        
    def ImportPoints(self, filename):
        points = np.genfromtxt(filename, dtype = 'float', delimiter='\t')
        
        self.px = np.asarray(points[:,0])
        self.py = np.asarray(points[:,1])
        self.pz = np.asarray(points[:,2])
        self.pstrike = np.asarray(points[:,3])
        self.pdip = np.asarray(points[:,4])
        self.prake = np.asarray(points[:,5])
        self.pm = np.asarray(points[:,6])
        
    def ImportBeta(self, filename):
        points = np.genfromtxt(filename, dtype = 'float', delimiter='\t')
        
        self.Beta = np.asarray(points[:,6])
        
    def CalculateBeta(self):
        regr = linear_model.LinearRegression(positive=True)
        regr.fit(self.X, self.y)
    
        self.Beta = np.asarray(regr.coef_)
        
    def CalculateBetaLasso(self):
        regr = linear_model.LassoLarsCV(positive=True)
        regr.fit(self.X, self.y)
        
        #print(regr.coef_)
    
        self.Beta = np.asarray(regr.coef_)
        
    def fnc2min(self, params, indices, data):
        coeffs = []
        for key, value in params.valuesdict().items():
            coeffs.append(value)
            
        model = []
        for i in indices:
            
            greensvals = self.X[i,:]

            model.append(np.sum(np.multiply(coeffs, greensvals)))

        return model - self.y

    def CalculateBetaLM(self, max=np.inf):
    
        indices = np.arange(0, len(self.y))
        
        params = Parameters()
        for j in range(self.X.shape[1]):
            name = "p%i" % j
            #val = 1.0e9
            val = np.random.uniform(1.0e6, 1.0e9)
            params.add(name, value=val, min=0, max = max)
            
        minner = Minimizer(self.fnc2min, params, fcn_args=(indices, self.y))
        result = minner.minimize()
        
        #print(fit_report(result))
        self.Beta = np.asarray([value for key, value in result.params.valuesdict().items()])

    def SaveModel(self, outfilename):
        
        with open(outfilename, 'w') as ofile:
            for px, py, pz, pstrike, pdip, prake, pm in zip(self.px, self.py, self.pz, self.pstrike, self.pdip, self.prake, self.Beta):
                ofile.write(str(px) + "\t" + str(py) + "\t" + str(pz) + "\t" + str(pstrike) + "\t" + str(pdip) + "\t" + str(prake) + "\t" + str(pm) + "\n")
                
    def CalculateBetaLMGuess(self, max=np.inf, guess=""):
    
        if not guess:
            self.CalculateBetaLM()
        else:
            points = np.genfromtxt(guess, dtype = 'float', delimiter='\t')
            self.Beta = np.asarray(points[:,6])

        indices = np.arange(0, len(self.y))
        
        params = Parameters()
        for j, val in zip(range(self.X.shape[1]), self.Beta):
            name = "p%i" % j
            params.add(name, value=val, min=0, max = max)
            
        minner = Minimizer(self.fnc2min, params, fcn_args=(indices, self.y))
        result = minner.minimize()
        
        #print(fit_report(result))
        self.Beta = np.asarray([value for key, value in result.params.valuesdict().items()])
    
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
        
    def PlotSmoothedData(self, plotid):
    
        self.dataid = plotid
    
        Zeros = np.zeros(len(self.datastore.data))
        tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(self.datastore.data[:,0], self.datastore.data[:,1], Zeros)
        
        self.ax[plotid].plot_trisurf(np.add(self.datastore.datasmooth[:,0], self.datastore.datasmooth[:,2]), self.datastore.datasmooth[:,1], self.datastore.datasmooth[:,3], triangles=tri.triangles, cmap = 'magma')
        self.vmin = self.datastore.datasmooth[:,3].min()
        self.vmax = self.datastore.datasmooth[:,3].max()
        self.ax[plotid].set_xlabel('x (km)')
        self.ax[plotid].set_ylabel('y (km)')
        self.ax[plotid].set_zlabel('z (m)')
        self.ax[plotid].set_title("Data")
        #self.ax[plotid].view_init(azim=45, elev=20)
        self.ax[plotid].view_init(azim=-90, elev=90)
        
    def PlotModel(self, plotid):
    
        Zeros = np.zeros(len(self.datastore.data))
        tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(self.datastore.data[:,0], self.datastore.data[:,1], Zeros)
        modeldata = [np.sum(np.multiply(self.datastore.Beta, x)) for x in self.datastore.X]
    
        self.ax[plotid].plot_trisurf(np.add(self.datastore.data[:,0], modeldata[:len(self.datastore.data)]), self.datastore.data[:,1], modeldata[len(self.datastore.data):], triangles=tri.triangles, cmap = 'magma', vmin = self.vmin, vmax = self.vmax)
        self.ax[plotid].set_xlabel('x (km)')
        self.ax[plotid].set_ylabel('y (km)')
        self.ax[plotid].set_zlabel('z (m)')
        self.ax[plotid].set_title("Model")
        #self.ax[plotid].view_init(azim=45, elev=20)
        self.ax[plotid].view_init(azim=-90, elev=90)
        
    def PlotResiduals(self, plotid):
    
        Zeros = np.zeros(len(self.datastore.data))
        tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(self.datastore.data[:,0], self.datastore.data[:,1], Zeros)
        modeldata = [np.sum(np.multiply(self.datastore.Beta, x)) for x in self.datastore.X]
    
        self.ax[plotid].plot_trisurf(np.add(self.datastore.data[:,0], modeldata[:len(self.datastore.data)]), self.datastore.data[:,1], modeldata[len(self.datastore.data):] - self.datastore.data[:,3], triangles=tri.triangles, cmap = 'magma', vmin = self.vmin, vmax = self.vmax)
        self.ax[plotid].set_xlabel('x (km)')
        self.ax[plotid].set_ylabel('y (km)')
        self.ax[plotid].set_zlabel('z (m)')
        self.ax[plotid].set_title("Residuals")
        #self.ax[plotid].view_init(azim=45, elev=20)
        self.ax[plotid].view_init(azim=-90, elev=90)
        
        
    def PlotModelPoints(self, plotid):
    
        px = self.datastore.px[np.where(self.datastore.Beta != 0)]
        py = self.datastore.py[np.where(self.datastore.Beta != 0)]
        pz = self.datastore.pz[np.where(self.datastore.Beta != 0)]
        m = self.datastore.Beta[np.where(self.datastore.Beta != 0)]
    
        #scatter = self.ax[plotid].scatter(self.datastore.px, self.datastore.py, self.datastore.pz, c = self.datastore.Beta, norm=colors.LogNorm(vmin=np.asarray(self.datastore.Beta).min(), vmax=np.asarray(self.datastore.Beta).max()), s = 100)
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
    
    store = DataStore()
    # store.ImportData("ridgecrest_displacements_inc_fitted_data.txt")
    # store.ImportGreens()
    # store.ImportPoints("ridgecrest_displacements_inc_model_points.txt")
    
    # store.ImportData("ridgecrest_binned_fitted_data.txt")
    # store.ImportGreens()
    # store.ImportPoints("ridgecrest_binned_model_points.txt")
    
    store.ImportData("iran_binned_fitted_data.txt")
    store.ImportGreens()
    store.ImportPoints("iran_binned_model_points.txt")
    
    start = time.time()
    #store.CalculateBeta()
    store.CalculateBetaLMGuess()
    store.SaveModel("optimizedmodel.txt")
    
    end = time.time()
    print(end - start)
    
    plotter = Plotter(store, 4)
    plotter.PlotData(0)
    #plotter.PlotSmoothedData(0)
    plotter.PlotModel(1)
    #plotter.PlotSmoothedData(1)
    plotter.PlotModelPoints(2)
    plotter.PlotResiduals(3)
    
    plt.show()