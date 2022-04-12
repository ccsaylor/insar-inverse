import numpy as np

class System:

    #generates the point sources in a model based on the parameters imported from sourceparams
    def __init__(self, sourceparams):
        
        sx = np.linspace(sourceparams.sxmin, sourceparams.sxmax, sourceparams.nsx)
        sy = np.linspace(sourceparams.symin, sourceparams.symax, sourceparams.nsy)
        sz = np.linspace(sourceparams.szmin, sourceparams.szmax, sourceparams.nsz)
        
        with open(sourceparams.filename[:-4] + "_model_points.txt", "w") as fout:
            for z in sz:
                for y in sy:
                    for x in sx:
                        fout.write(str(x) + '\t')
                        fout.write(str(y) + '\t')
                        fout.write(str(z) + '\t')
                        fout.write(str(sourceparams.strike) + '\t')
                        fout.write(str(sourceparams.dip) + '\t')
                        fout.write(str(sourceparams.rake) + '\n')

    
    