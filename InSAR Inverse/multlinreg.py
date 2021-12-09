from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import pandas as pd
import matplotlib.tri as mtri
from sklearn import linear_model
from lmfit import Minimizer, Parameters, fit_report
import time
    
class DataStore:
    def __init__(self):
        self.test = []
    
    def ImportData(self, filename, xmin = -np.inf, xmax = np.inf, ymin = -np.inf, ymax = np.inf):
        colnames = ["x", "y", "z"]
        df = pd.read_table(filename, names=colnames).dropna()
        self.y = df[(df.x > xmin) & (df.x < xmax) & (df.y > ymin) & (df.y < ymax)]
        self.ynumpy = self.y.to_numpy()
        
    def ImportGreens(self, filename):
        matrix = np.genfromtxt(filename, dtype = 'float', delimiter='\t')
        
        colnames = []
        for i, val in enumerate(matrix[0]):
            colnames.append("Source %i" % i)

        self.X = pd.read_table(filename, names=colnames)
        self.Xnumpy = self.X.to_numpy()
        
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
        regr.fit(self.X, self.y["z"])
    
        self.Beta = np.asarray(regr.coef_)
        
    # def CalculateBetaRidge(self): #doesn't work super well
        # regr = linear_model.Ridge(alpha=0.5)
        # regr.fit(self.X, self.y["z"])
    
        # self.Beta = np.asarray(regr.coef_)
        
    def CalculateBetaLasso(self):
        regr = linear_model.LassoLarsCV(positive=True)
        regr.fit(self.X, self.y["z"])
        
        print(regr.coef_)
    
        self.Beta = np.asarray(regr.coef_)
        
    # def CalculateBetaElasticNet(self): #not sure how to make this work
        # regr = linear_model.ElasticNetCV(positive=True)
        # regr.fit(self.X, self.y["z"])
    
        # self.Beta = np.asarray(regr.coef_)
        
    def fnc2min(self, params, indices, data):
        coeffs = []
        for key, value in params.valuesdict().items():
            coeffs.append(value)
            
        model = []
        for i in indices:
            
            greensvals = self.Xnumpy[i,:]

            model.append(np.sum(np.multiply(coeffs, greensvals)))

        return model - data[:,2]

    def CalculateBetaLM(self, max=np.inf):
    
        indices = np.arange(0, len(self.y))
        
        params = Parameters()
        for j in range(self.X.shape[1]):
            name = "p%i" % j
            params.add(name, value=1.0e9, min=0, max = max)
            
        minner = Minimizer(self.fnc2min, params, fcn_args=(indices, self.ynumpy))
        result = minner.minimize()
        
        #print(fit_report(result))
        self.Beta = [value for key, value in result.params.valuesdict().items()]

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
            
        minner = Minimizer(self.fnc2min, params, fcn_args=(indices, self.ynumpy))
        result = minner.minimize()
        
        #print(fit_report(result))
        self.Beta = [value for key, value in result.params.valuesdict().items()]
    
# class Plotter:
    
    # def __init__(self, store, num_plots):
        # self.datastore = store
        # self.fig, self.ax = plt.subplots(1, num_plots, facecolor='w', edgecolor='k', subplot_kw=dict(projection='3d'))
        # self.dataid = -1

    # def PlotData(self, plotid):
    
        # self.dataid = plotid
    
        # Zeros = np.zeros(len(self.datastore.y["z"]))
        # tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(self.datastore.y["x"], self.datastore.y["y"], Zeros)
        
        # self.vmin = self.datastore.y["z"].min()
        # self.vmax = self.datastore.y["z"].max()
        
        # self.ax[plotid].plot_trisurf(self.datastore.y["x"], self.datastore.y["y"], self.datastore.y["z"], triangles=tri.triangles, cmap = 'magma')
        # self.ax[plotid].set_xlabel('x (km)')
        # self.ax[plotid].set_ylabel('y (km)')
        # self.ax[plotid].set_zlabel('z (m)')
        # self.ax[plotid].set_title("Data")
        # #self.ax[plotid].view_init(azim=45, elev=20)
        # self.ax[plotid].view_init(azim=-90, elev=90)
        # self.ax[plotid].set_box_aspect((np.ptp(self.datastore.y["x"]), np.ptp(self.datastore.y["y"]), np.ptp(self.datastore.y["z"])))
        
    # def PlotModel(self, plotid):
    
        # Zeros = np.zeros(len(self.datastore.y["z"]))
        # tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(self.datastore.y["x"], self.datastore.y["y"], Zeros)
        # modeldata = [np.sum(np.multiply(self.datastore.Beta, x)) for index, x in self.datastore.X.iterrows()]
    
        # self.ax[plotid].plot_trisurf(self.datastore.y["x"], self.datastore.y["y"], modeldata, triangles=tri.triangles, cmap = 'magma', vmin = self.vmin, vmax = self.vmax)
        # self.ax[plotid].set_xlabel('x (km)')
        # self.ax[plotid].set_ylabel('y (km)')
        # self.ax[plotid].set_zlabel('z (m)')
        # self.ax[plotid].set_title("Model")
        # #self.ax[plotid].view_init(azim=45, elev=20)
        # self.ax[plotid].view_init(azim=-90, elev=90)
        
    # def PlotResiduals(self, plotid):
    
        # Zeros = np.zeros(len(self.datastore.y["z"]))
        # tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(self.datastore.y["x"], self.datastore.y["y"], Zeros)
        # modeldata = [np.sum(np.multiply(self.datastore.Beta, x)) for index, x in self.datastore.X.iterrows()]
    
        # self.ax[plotid].plot_trisurf(self.datastore.y["x"], self.datastore.y["y"], modeldata - self.datastore.y["z"], triangles=tri.triangles, cmap = 'magma', vmin = self.vmin, vmax = self.vmax)
        # self.ax[plotid].set_xlabel('x (km)')
        # self.ax[plotid].set_ylabel('y (km)')
        # self.ax[plotid].set_zlabel('z (m)')
        # self.ax[plotid].set_title("Residuals")
        # #self.ax[plotid].view_init(azim=45, elev=20)
        # self.ax[plotid].view_init(azim=-90, elev=90)
        
    # def PlotModelPoints(self, plotid):
    
        # px = self.datastore.px[np.where(self.datastore.Beta != 0)]
        # py = self.datastore.py[np.where(self.datastore.Beta != 0)]
        # pz = self.datastore.pz[np.where(self.datastore.Beta != 0)]
        # m = self.datastore.Beta[np.where(self.datastore.Beta != 0)]
    
        # #scatter = self.ax[plotid].scatter(self.datastore.px, self.datastore.py, self.datastore.pz, c = self.datastore.Beta, norm=colors.LogNorm(vmin=np.asarray(self.datastore.Beta).min(), vmax=np.asarray(self.datastore.Beta).max()), s = 100)
        # scatter = self.ax[plotid].scatter(px, py, pz, c = m, s = 20)

        # self.ax[plotid].set_xlabel('x (km)')
        # self.ax[plotid].set_ylabel('y (km)')
        # self.ax[plotid].set_zlabel('z (m)')
        # #self.ax[plotid].view_init(azim=45, elev=20)
        # self.ax[plotid].view_init(azim=-90, elev=90)
        
        # if self.dataid != -1:
            # self.ax[plotid].set_xlim(self.ax[self.dataid].get_xlim())
            # self.ax[plotid].set_ylim(self.ax[self.dataid].get_ylim())
        # #plt.colorbar(scatter)
        
class Plotter:
    
    def __init__(self, store, nrows, ncols):
        self.datastore = store
        self.nrows = nrows
        self.ncols = ncols
        self.fig, self.ax = plt.subplots(nrows, ncols, facecolor='w', edgecolor='k', subplot_kw=dict(projection='3d'))
        self.dataid = -1

    def PlotData(self, plotid):
        self.dataid = plotid
        
        i = plotid // self.nrows
        j = plotid % self.ncols
    
        Zeros = np.zeros(len(self.datastore.y["z"]))
        tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(self.datastore.y["x"], self.datastore.y["y"], Zeros)
        
        self.vmin = self.datastore.y["z"].min()
        self.vmax = self.datastore.y["z"].max()
        
        self.ax[i, j].plot_trisurf(self.datastore.y["x"], self.datastore.y["y"], self.datastore.y["z"], triangles=tri.triangles, cmap = 'magma')
        self.ax[i, j].set_xlabel('x (km)')
        self.ax[i, j].set_ylabel('y (km)')
        self.ax[i, j].set_zlabel('z (m)')
        self.ax[i, j].set_title("Data")
        #self.ax[i, j].view_init(azim=45, elev=20)
        self.ax[i, j].view_init(azim=-90, elev=90)
        self.ax[i, j].set_box_aspect((np.ptp(self.datastore.y["x"]), np.ptp(self.datastore.y["y"]), np.ptp(self.datastore.y["z"])))
        
    def PlotModel(self, plotid):
        i = plotid // self.nrows
        j = plotid % self.ncols
    
        Zeros = np.zeros(len(self.datastore.y["z"]))
        tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(self.datastore.y["x"], self.datastore.y["y"], Zeros)
        modeldata = [np.sum(np.multiply(self.datastore.Beta, x)) for index, x in self.datastore.X.iterrows()]
    
        self.ax[i, j].plot_trisurf(self.datastore.y["x"], self.datastore.y["y"], modeldata, triangles=tri.triangles, cmap = 'magma', vmin = self.vmin, vmax = self.vmax)
        self.ax[i, j].set_xlabel('x (km)')
        self.ax[i, j].set_ylabel('y (km)')
        self.ax[i, j].set_zlabel('z (m)')
        self.ax[i, j].set_title("Model")
        #self.ax[i, j].view_init(azim=45, elev=20)
        self.ax[i, j].view_init(azim=-90, elev=90)
        self.ax[i, j].set_box_aspect((np.ptp(self.datastore.y["x"]), np.ptp(self.datastore.y["y"]), np.ptp(modeldata)))
        
    def PlotResiduals(self, plotid):
        i = plotid // self.nrows
        j = plotid % self.ncols
    
        Zeros = np.zeros(len(self.datastore.y["z"]))
        tri, args, kwargs = mtri.Triangulation.get_from_args_and_kwargs(self.datastore.y["x"], self.datastore.y["y"], Zeros)
        modeldata = [np.sum(np.multiply(self.datastore.Beta, x)) for index, x in self.datastore.X.iterrows()]
    
        self.ax[i, j].plot_trisurf(self.datastore.y["x"], self.datastore.y["y"], modeldata - self.datastore.y["z"], triangles=tri.triangles, cmap = 'magma', vmin = self.vmin, vmax = self.vmax)
        self.ax[i, j].set_xlabel('x (km)')
        self.ax[i, j].set_ylabel('y (km)')
        self.ax[i, j].set_zlabel('z (m)')
        self.ax[i, j].set_title("Residuals")
        #self.ax[i, j].view_init(azim=45, elev=20)
        self.ax[i, j].view_init(azim=-90, elev=90)
        self.ax[i, j].set_box_aspect((np.ptp(self.datastore.y["x"]), np.ptp(self.datastore.y["y"]), np.ptp(modeldata - self.datastore.y["z"])))
        
    def PlotModelPoints(self, plotid):
        i = plotid // self.nrows
        j = plotid % self.ncols
    
        px = self.datastore.px[np.where(self.datastore.Beta != 0)]
        py = self.datastore.py[np.where(self.datastore.Beta != 0)]
        pz = self.datastore.pz[np.where(self.datastore.Beta != 0)]
        m = self.datastore.Beta[np.where(self.datastore.Beta != 0)]
    
        #scatter = self.ax[i, j].scatter(self.datastore.px, self.datastore.py, self.datastore.pz, c = self.datastore.Beta, norm=colors.LogNorm(vmin=np.asarray(self.datastore.Beta).min(), vmax=np.asarray(self.datastore.Beta).max()), s = 100)
        scatter = self.ax[i, j].scatter(px, py, pz, c = m, s = 20)

        self.ax[i, j].set_xlabel('x (km)')
        self.ax[i, j].set_ylabel('y (km)')
        self.ax[i, j].set_zlabel('z (m)')
        #self.ax[i, j].view_init(azim=45, elev=20)
        self.ax[i, j].view_init(azim=-90, elev=90)
        self.ax[i, j].set_box_aspect((np.ptp(self.datastore.y["x"]), np.ptp(self.datastore.y["y"]), np.ptp(self.datastore.y["z"])))
        
        if self.dataid != -1:
            id = self.dataid // self.nrows
            jd = self.dataid % self.ncols
            self.ax[i, j].set_xlim(self.ax[id, jd].get_xlim())
            self.ax[i, j].set_ylim(self.ax[id, jd].get_ylim())
        #plt.colorbar(scatter)
        
if __name__ == "__main__":
    
    store = DataStore()
    store.ImportData("nepal_binned_75x75_converted_fitted_data.txt")
    store.ImportGreens("greens_z.txt")
    store.ImportPoints("nepal_binned_75x75_converted_model_points.txt")
    
    start = time.time()
    #store.ImportBeta("39x27optimizedfull.txt")
    store.CalculateBeta()
    #store.CalculateBetaLMGuess(guess="26x18guess.txt")
    store.SaveModel("optimizedmodel.txt")
    
    end = time.time()
    print(end - start)
    
    plotter = Plotter(store, 2, 2)
    plotter.PlotData(0)
    plotter.PlotModel(1)
    plotter.PlotModelPoints(2)
    plotter.PlotResiduals(3)
    
    plt.show()
    