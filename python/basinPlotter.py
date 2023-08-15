import basinFunctions
from config import Config 
import numpy as np
import glob
import os
import imageio

config = Config()
phaseShifts = np.linspace(0,np.pi,10)
for ph in phaseShifts:
    config.phaseShift = ph
    
    

if __name__ == "__main__":
    photoOrderList = np.linspace(0,359,360, dtype=int)
    fileList = []
    print(photoOrderList)
    #files = glob.glob(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\alongSQRT2RadCircle\*.csv")
    #for photo in photoOrderList[262:]:
    #    print(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\alongSQRT2RadCircle\*_"+str(photo)+r"_["+r"*.csv")
    #    file = glob.glob(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\alongSQRT2RadCircle\*_"+str(photo)+r"_["+r"*.csv")[0]
    #    simData = basinFunctions.readCSV(file)
    #    basinFunctions.plot_newton_fractal(simData, config, file)
    
    basinFunctions.gifBasins(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\alongSQRT2RadCircle\\")


