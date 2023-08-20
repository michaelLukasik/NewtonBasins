import numpy as np 
import matplotlib.pyplot as plt
import glob


piDomain = "[-3.141590_3.141590_-3.141590_3.141590_]"
pi2Domain= "[-1.570795_1.570795_-1.570795_1.570795_]"
pi4Domain = "[-0.785397_0.785397_-0.785397_0.785397_]"
pi16Domain = "[-0.196349_0.196349_-0.196349_0.196349_]"
pi10Domain = "[-0.314159_0.314159_-0.314159_0.314159_]"
t100piDomain = "[-314.159000_314.159000_-314.159000_314.159000_]" 
t10piDomain = "[-31.415900_31.415900_-31.415900_31.415900_]"
t4piDomain = "[-12.566360_12.566360_-12.566360_12.566360_]"
oneDomain = "[-1.000000_1.000000_-1.000000_1.000000_]"
tenDomain = "[-10.000000_10.000000_-10.000000_10.000000_]"


##Bessel
offset = "[3.141590,0.000000]"
class Config: 
    
    colors =  ['b', 'r', 'g', 'y']
    n = 50 ## Number of divisions on the screen (n x n total pixels)
    domain = [-np.pi, np.pi, -np.pi, np.pi]  
    tol = 1e-10
    cmap = "hsv"
    its = "100"
    #file = r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\CustomDriver1\\" +tenDomain+r"_"+its+"_"+str(n)+".csv"
    file = r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\\variableOffset\\VaryingReals\\" +t4piDomain+r"_"+its+"_"+str(n)+"_"+offset+".csv"
    phaseShift = 0.
    
    
def makeBasin(m, config, file):
    print("Making basin")
    plt.figure(figsize = (10,10))
    plt.imshow(m, cmap=config.cmap, origin='lower')
    plt.axis('off')
    plt.savefig(file[:-4] + "_"+ config.cmap + "_" + str(config.tol) + ".png",bbox_inches='tight', transparent="True", pad_inches=0)
    #plt.show() 