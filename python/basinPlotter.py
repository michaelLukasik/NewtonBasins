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
    photoOrderList = np.linspace(0,0,1, dtype=int)
    fileList = []
    print(photoOrderList)
    for photo in photoOrderList:
        file = glob.glob(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\CoshSinc\variableOffset\\Halley\\*_30_3000_*"+str(photo)+r"_["+r"*.csv")[0]
        simData = basinFunctions.readCSV(file)
        basinFunctions.plot_newton_fractal(simData, config, file)    
    
    
    ### Path 1 
    #photo0List = np.linspace(0,99,100, dtype=int) ## Up the quad
    #photo1List = np.linspace(99,30,70, dtype=int) ## Down the quad
    #photo2List = np.linspace(0,179,180, dtype=int) ## Around Circle
    #photo3List = np.linspace(0,179,180, dtype=int) ## Down Spiral
    
    #for photo in photo0List:
    #     path0File = glob.glob(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\VaryingAlongQuadExtended\\" + r"*_"+str(photo)+r"_["+r"*.png")
    #     fileList.append(path0File)
    #for photo in photo1List:
    #    path1File = glob.glob(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\VaryingAlongQuadExtended\\" + r"*_"+str(photo)+r"_["+r"*.png")
    #    fileList.append(path1File)
    #for photo in photo2List:
    #    path2File = glob.glob(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\VaryingAlongUnitCircle\\" + r"*_"+str(photo)+r"_["+r"*.png")
    #    fileList.append(path2File)
    #for photo in photo3List:
    #    path3File = glob.glob(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\AlongSpiralInwards\\" + r"*_"+str(photo)+r"_["+r"*.png")
    #    fileList.append(path3File)

    #print(fileList)
    
    ### Path 2 Big Spiral Inwards
    
    #photoList = np.linspace(0,359,360,dtype=int)
    #for photo in photoList:
    #     path0File = glob.glob(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\SpiralInwards4Cycle2Radius360Images\\" + r"*_"+str(photo)+r"_["+r"*.png")
    #     fileList.append(path0File)
    
    ### Path 3 LJ -- 
    
    photoList = np.linspace(0,89,90,dtype=int)
    for photo in photoList:
        fileName = "LissaJous_A2B3D3pi4"
        path0File = glob.glob(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\\"+fileName+"\\" + r"*_"+str(photo)+r"_["+r"*.png")
        fileList.append(path0File)
    
    basinFunctions.gifBasinsFileList(fileList, fileName)

