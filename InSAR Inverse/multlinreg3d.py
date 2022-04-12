from mpl_toolkits.mplot3d import Axes3D
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import matplotlib.tri as mtri
from sklearn import linear_model
from lmfit import Minimizer, Parameters, fit_report
import time
from pointokada import PointOkada
from sourceparams import SourceParams
from system import System
from scipy.interpolate import griddata
    
class DataStore:
    def __init__(self, usedx, usedy, usedz):
        self.usedx = usedx
        self.usedy = usedy
        self.usedz = usedz
    
    def ImportData(self, filename):
        self.data = np.genfromtxt(filename, dtype = 'float', delimiter='\t')
        ncols = self.data.shape[1]
        self.y = []
        for i in range(2, ncols):
            self.y = np.concatenate((self.y, self.data[:,i]), axis=None)
            
        self.xmin = self.data[:,0].min()
        self.xmax = self.data[:,0].max()
        self.ymin = self.data[:,1].min()
        self.ymax = self.data[:,1].max()
        
    def GeneratePtsGreens(self, paramfile):
        params = SourceParams(paramfile, self)
        sys = System(params)
    
        points = np.genfromtxt(paramfile[:-4] + "_model_points.txt", dtype = 'float', delimiter='\t')
        
        self.px = np.asarray(points[:,0])
        self.py = np.asarray(points[:,1])
        self.pz = np.asarray(points[:,2])
        self.pstrike = np.asarray(points[:,3])
        self.pdip = np.asarray(points[:,4])
        self.prake = np.asarray(points[:,5])
        
        self.X = np.zeros(1)
        
        po = PointOkada()
        if self.usedx:
            with open("greens_x.txt", "w") as fout:
                for x, y in self.data[:,:2]:
                    for i, (px, py, pz, strike, dip, rake) in enumerate(zip(self.px, self.py, self.pz, self.pstrike, self.pdip, self.prake)):
                        dx = po.CalcGreensX(x, y, px, py, np.abs(pz), strike, dip, rake)
                        if i == len(self.px) - 1:
                            fout.write(str(dx))
                        else:
                            fout.write(str(dx) + '\t')
                    fout.write('\n')
            ew = np.genfromtxt("greens_x.txt", dtype = 'float', delimiter='\t')
            
            if not self.X.all():
                self.X = ew
            else:
                self.X = np.concatenate((self.X, ew), axis=0)
        
        if self.usedy:
            with open("greens_y.txt", "w") as fout:
                for x, y in self.data[:,:2]:
                    for i, (px, py, pz, strike, dip, rake) in enumerate(zip(self.px, self.py, self.pz, self.pstrike, self.pdip, self.prake)):
                        dy = po.CalcGreensY(x, y, px, py, np.abs(pz), strike, dip, rake)
                        if i == len(self.px) - 1:
                            fout.write(str(dy))
                        else:
                            fout.write(str(dy) + '\t')
                    fout.write('\n')
            ns = np.genfromtxt("greens_y.txt", dtype = 'float', delimiter='\t')
            
            if not self.X.all():
                self.X = ns
            else:
                self.X = np.concatenate((self.X, ns), axis=0)
                    
        if self.usedz:
            with open("greens_z.txt", "w") as fout:
                for x, y in self.data[:,:2]:
                    for i, (px, py, pz, strike, dip, rake) in enumerate(zip(self.px, self.py, self.pz, self.pstrike, self.pdip, self.prake)):
                        dz = po.CalcGreensZ(x, y, px, py, np.abs(pz), strike, dip, rake)
                        if i == len(self.px) - 1:
                            fout.write(str(dz))
                        else:
                            fout.write(str(dz) + '\t')
                    fout.write('\n')
            ud = np.genfromtxt("greens_z.txt", dtype = 'float', delimiter='\t')
            
            if not self.X.all():
                self.X = ud
            else:
                self.X = np.concatenate((self.X, ud), axis=0)
            
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
    
    def __init__(self, store):
        self.datastore = store
        
    def PlotComparison(self):

        pts = self.datastore.data[:,:2]
        data = np.copy(self.datastore.y)
        modeldata = [np.sum(np.multiply(self.datastore.Beta, x)) for x in self.datastore.X]
        
        xmin = self.datastore.xmin
        xmax = self.datastore.xmax
        ymin = self.datastore.ymin
        ymax = self.datastore.ymax
        num = 30j
        
        x, y = np.mgrid[xmin:xmax:num, ymin:ymax:num]
        
        vmin = data.min()
        vmax = data.max()
        
        print(vmin, vmax)
        
        if self.datastore.usedx:
            figx, axx = plt.subplots(1, 2)
            nrows = pts.shape[0]
            datadx = data[:nrows]
            modeldx = modeldata[:nrows]
            resampleddata = griddata((pts[:,0], pts[:,1]), datadx, (x, y))
            resampledmodel = griddata((pts[:,0], pts[:,1]), modeldx, (x, y))
            
            axx[0].imshow(resampleddata.T, origin='lower', cmap='coolwarm', extent=(xmin, xmax, ymin, ymax), vmin = vmin, vmax = vmax)
            axx[1].imshow(resampledmodel.T, origin='lower', cmap='coolwarm', extent=(xmin, xmax, ymin, ymax), vmin = vmin, vmax = vmax)
            axx[0].set_title("Data dx")
            axx[1].set_title("Model dx")
            
            data = data[nrows:]
            modeldata = modeldata[nrows:]
            
        if self.datastore.usedy:
            figy, axy = plt.subplots(1, 2)
            nrows = pts.shape[0]
            datady = data[:nrows]
            modeldy = modeldata[:nrows]
            resampleddata = griddata((pts[:,0], pts[:,1]), datady, (x, y))
            resampledmodel = griddata((pts[:,0], pts[:,1]), modeldy, (x, y))
            
            axy[0].imshow(resampleddata.T, origin='lower', cmap='coolwarm', extent=(xmin, xmax, ymin, ymax), vmin = vmin, vmax = vmax)
            axy[1].imshow(resampledmodel.T, origin='lower', cmap='coolwarm', extent=(xmin, xmax, ymin, ymax), vmin = vmin, vmax = vmax)
            axy[0].set_title("Data dy")
            axy[1].set_title("Model dy")
            
            data = data[nrows:]
            modeldata = modeldata[nrows:]
            
        if self.datastore.usedz:
            figz, axz = plt.subplots(1, 2)
            nrows = pts.shape[0]
            datadz = data[:nrows]
            modeldz = modeldata[:nrows]
            resampleddata = griddata((pts[:,0], pts[:,1]), datadz, (x, y))
            resampledmodel = griddata((pts[:,0], pts[:,1]), modeldz, (x, y))
            
            axz[0].imshow(resampleddata.T, origin='lower', cmap='coolwarm', extent=(xmin, xmax, ymin, ymax), vmin = vmin, vmax = vmax)
            axz[1].imshow(resampledmodel.T, origin='lower', cmap='coolwarm', extent=(xmin, xmax, ymin, ymax), vmin = vmin, vmax = vmax)
            axz[0].set_title("Data dz")
            axz[1].set_title("Model dz")
        
    def PlotModelPointsSolo(self):
    
        fig, ax = plt.subplots(1)
        
        px = self.datastore.px[np.where(self.datastore.Beta != 0)]
        py = self.datastore.py[np.where(self.datastore.Beta != 0)]
        m = self.datastore.Beta[np.where(self.datastore.Beta != 0)]
    
        ax.scatter(px, py, c = m, cmap = 'bone_r')
        ax.set_aspect("equal")
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')


        
if __name__ == "__main__":
    
    store = DataStore(True, False, True)
    
    start = time.time()
    
    store.ImportData("iran_binned_fitted_data.txt")
    store.GeneratePtsGreens("iranparams.ini")
    
    # store.ImportData("nepal_binned_75x75_converted_fitted_data.txt")
    # store.GeneratePtsGreens("nepalparams.ini")
    
    # store.ImportData("ridgecrest_binned_fitted_data.txt")
    # store.GeneratePtsGreens("ridgecrestbinnedparams.ini")
    
    # store.ImportData("ridgecrest_combined.txt")
    # store.GeneratePtsGreens("ridgecrestcombinedparams.ini")
    
    # store.ImportData("ridgecrest_100x100_kriged_thresh.txt")
    # store.GeneratePtsGreens("ridgecrestkrigedparams.ini")
    
    #store.CalculateBeta()
    store.CalculateBetaLMGuess()
    
    #store.SaveModel("optimizedmodel.txt")
    
    #store.ImportBeta("optimizedmodel.txt")
    
    end = time.time()
    print(end - start)
    
    plotter = Plotter(store)
    plotter.PlotComparison()
    plotter.PlotModelPointsSolo()
    
    plt.show()
