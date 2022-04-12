import configparser

class SourceParams:
    def __init__(self, paramfile, datastore):
        config = configparser.ConfigParser()
        config.read(paramfile)
        
        self.filename = paramfile

        self.xmin = config.getfloat("Data parameters", "Minimum data x value")
        self.xmax = config.getfloat("Data parameters", "Maximum data x value")
        self.ymin = config.getfloat("Data parameters", "Minimum data y value")
        self.ymax = config.getfloat("Data parameters", "Maximum data y value")
        
        self.strike = config.getfloat("Source parameters", "Strike angle")
        self.dip = config.getfloat("Source parameters", "Dip angle")
        self.rake = config.getfloat("Source parameters", "Rake angle")
        self.mmin = config.getfloat("Source parameters", "Minimum seismic moment")
        self.mmax = config.getfloat("Source parameters", "Maximum seismic moment")
        
        self.nsx = config.getint("Source parameters", "Sources in x direction")
        self.nsy = config.getint("Source parameters", "Sources in y direction")
        self.nsz = config.getint("Source parameters", "Sources in z direction")
        
        self.sxmin = config.getfloat("Source parameters", "Minimum source x value")
        self.sxmax = config.getfloat("Source parameters", "Maximum source x value")
        self.symin = config.getfloat("Source parameters", "Minimum source y value")
        self.symax = config.getfloat("Source parameters", "Maximum source y value")
        self.szmin = config.getfloat("Source parameters", "Minimum source z value")
        self.szmax = config.getfloat("Source parameters", "Maximum source z value")
        
        if self.xmin == -1 or self.xmax == -1 or self.ymin == -1 or self.ymax == -1:
            print("Data limits not set; using all data points")
            self.xmin = datastore.data[:,0].min()
            self.xmax = datastore.data[:,0].max()
            self.ymin = datastore.data[:,1].min()
            self.ymax = datastore.data[:,1].max()

        if self.sxmin == -1 or self.sxmax == -1 or self.symin == -1 or self.symax == -1:
            print("Source limits not set; setting to data limits")
            self.sxmin = self.xmin
            self.sxmax = self.xmax
            self.symin = self.ymin
            self.symax = self.ymax

    
    