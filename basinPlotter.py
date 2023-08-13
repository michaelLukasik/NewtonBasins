
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import ListedColormap
import math


def readCSV(file):
    df = pd.read_csv(file, header = None , dtype = np.float64, na_values = "-nan(ind)")
    df = df.fillna(0) 
    sim_wfRe =  df[0]
    sim_wfIm =  df[1]
    sim_xpos = df[2] 
    sim_ypos = df[3]
    return [sim_wfRe, sim_wfIm, sim_xpos, sim_ypos]


def get_root_index(roots, r):
    try:
        return np.where(np.isclose(roots, r, atol=1e-5))[0][0]
    except IndexError:
        roots.append(r)
        return len(roots) - 1

def plot_newton_fractal(simData, n, domain, TOL = 1e-5):
    
    roots = []
    m = np.zeros((n,n))
    z = simData[0] + 1j*simData[1]
    for iz0,z0 in enumerate(z):
        ir = get_root_index(roots, z0)
        m[math.floor(iz0/n), math.floor(iz0 % n)] = ir  
    nroots = len(roots)
    if nroots > len(colors):
        # Use a "continuous" colormap if there are too many roots.
        cmap = 'magma'
    else:
        # Use a list of colors for the colormap: one for each root.
        cmap = ListedColormap(colors[:nroots])
    m= m.reshape((n,n))
    print(roots)
    plt.imshow(m, cmap=cmap, origin='lower')
    plt.axis('off')
    plt.show()   

file = r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Sinh\1.csv"
colors = ['b', 'r', 'g', 'y']

simData = readCSV(file)
plot_newton_fractal(simData, n=1000, domain=(-np.pi, np.pi, -np.pi, np.pi), TOL = 1e-5)
