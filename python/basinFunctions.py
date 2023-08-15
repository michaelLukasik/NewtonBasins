import pandas as pd
import numpy as np
from PIL import Image
import math
import glob
from matplotlib.colors import ListedColormap
import basinFunctions
import matplotlib.pyplot as plt
import imageio


def readCSV(file):
    df = pd.read_csv(file, header = None , dtype = np.float64, na_values = "-nan(ind)")
    df = df.fillna(0) 
    sim_wfRe =  df[0]
    sim_wfIm =  df[1]
    sim_xpos = df[2] 
    sim_ypos = df[3]
    iteration = df[4]
    return [sim_wfRe, sim_wfIm, sim_xpos, sim_ypos, iteration]

def get_iteration_index(iterations, r, config):
    iterations.append(r)
    return iterations

def get_root_index(roots, r, config):
    try:
        return np.where(np.isclose(roots, r, atol=config.tol))[0][0]
    except IndexError:
        roots.append(r)
        return len(roots) - 1
    
def plot_newton_fractal(simData, config, file):  ## Adapted from the SciPy Docs on Newton's Method
    
    iterations = []
    roots = []
    m = np.zeros((config.n,config.n))
    z = simData[0] + 1j*simData[1]
    for iz0,z0 in enumerate(z):
        i_iter = basinFunctions.get_iteration_index(iterations, simData[4], config)
        m[math.floor(iz0/config.n), math.floor(iz0 % config.n)] = i_iter[0][iz0]
    nroots = len(roots)
    if nroots > len(config.colors):
        # Use a "continuous" colormap if there are too many roots.
        cmap = config.cmap
    else:
        # Use a list of colors for the colormap: one for each root.
        cmap = ListedColormap(config.colors[:nroots])
    m= m.reshape((config.n,config.n))
    makeBasin(m, config, file)
    
def makeBasin(m, config, file):
    print("Making basin")
    plt.figure(figsize = (10,10))
    plt.imshow(m, cmap=config.cmap, origin='lower')
    plt.axis('off')
    plt.savefig(file[:-4] + "_"+ config.cmap + "_" + str(config.tol) + ".png",bbox_inches='tight', transparent="True", pad_inches=0)
    #plt.show() 
   
def gifBasins(folder):
    photoOrderList = np.linspace(0,49,50, dtype=int)
    frames = []
    for photo in photoOrderList:
        print(folder + "*_"+str(photo)+r"_["+r"*.png")
        imageFile = glob.glob(folder + "*_"+str(photo)+r"_["+r"*.png")[0]
        print(photo, imageFile)
        new_frame = Image.open(imageFile)
        frames.append(new_frame)
    imageio.mimsave(folder + 'GIF.gif', frames, loop = 0, duration = 1000 * 1/120)
    #frames[0].save(folder + 'GIF.gif', format='GIF',
    #           append_images=frames[1:],
    #           save_all=True,
    #           fps=120, loop=True)
    
    
def gifBasinsFileList(List):
    frames = []
    for photo in List:
        print(photo)
        new_frame = Image.open(photo[0])
        frames.append(new_frame)
    imageio.mimsave(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\Path1GIF.gif", frames, loop = 0, duration = 20)