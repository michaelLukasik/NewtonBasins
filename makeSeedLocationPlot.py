import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import imageio
from PIL import Image, ImageSequence, ImageDraw, ImageFont
import sympy

besselFormula = r"$p(z) = \frac{sin(z)}{z} + c$"
#sympy.preview(besselFormula, viewer='file', filename='besselFormula.png', euler=False, dvioptions=["-T", "tight", "-z", "0", "--truecolor", "-D 600", "-bg", "Transparent"])


writergif = animation.PillowWriter(fps=50)
bigfnt = ImageFont.truetype(r"C:\Windows\Fonts\cambria.ttc", 40)
fnt = ImageFont.truetype(r"C:\Windows\Fonts\cambria.ttc", 25)

zsymb = sympy.symbols('z')
expr = sympy.sin(zsymb) /zsymb 
lat = sympy.latex(expr) 



def LJLocation(a,b,d,index): 
    return a*np.cos( 2*np.pi *4* (index/360) + d) , b*np.sin( 2*np.pi *4* (index/360))

def animate(i):
    x,y = pathLocation(i)
    seedLocation.set_data([x[i],y[i]])
    seedLocation.set_label("c = ("+str(round(x[i],3))+","+str(round(y[i],3))+")")
    #seedLocation.set_alpha(1-(i*3e-3))
    plt.legend(loc = 'lower left')
    return seedLocation

def pathLocation (i):
    return np.array([fullPathX,fullPathY])

def examineGifs():
    gif1 = imageio.get_reader('Path2.gif')
    gif2 = imageio.get_reader(r"C:\Users\Michael\Documents\Programming\NewtonBasins\Path2.gif")
    print(gif1.get_length())
    print(gif2.get_length())
    
    
def embedGifs(pathString):
    equation = Image.open(r"Share\\BesselEqn.png")
    path   =  Image.open(r"paths\\"+pathString+r"Path.gif")
    pattern =  Image.open(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\\"+str(pathString)+"\\"+str(pathString)+".gif")
    totalH = max(path.size[0] ,pattern.size[0])
    totalW = max(path.size[1] ,equation.size[1]) + pattern.size[1]
    frames = []
    for i in range(391):
        print("path n frames " ,path.n_frames)
        print(i , pattern.n_frames)
        path.seek(i % path.n_frames)
        pattern.seek(i % pattern.n_frames)
        frame = Image.new('RGB', (totalW, totalH), (255, 255, 255))
        frame.paste(path, (0,0), mask=path.convert("RGBA"))
        frame.paste(pattern, (max(path.size[0] ,equation.size[0]),0), mask=pattern.convert("RGBA"))
        frame.paste(equation, (60,550) , mask=equation.convert("RGBA"))
        frames.append(frame)
    frames[0].save(str(pathString)+r".gif", save_all=True, duration=20, loop=0, append_images=frames[1:])


  #####################################################################  
    
fullPathX = []
fullPathY = []

'''  Path 1
unitCircleX = [np.cos(x) for x in np.linspace(np.pi/4.,2*np.pi+np.pi/4.,180)] 
unitCircleY = [np.sin(y) for y in np.linspace(np.pi/4.,2*np.pi+np.pi/4.,180)] 

parabolaX = [(x / 35) for x in np.linspace(0,99,100)]
parabolaY = [(y / 35)**2 for y in np.linspace(0,99,100)]

parabolaRetX = [(x / 35) for x in np.linspace(99,30,70)]
parabolaRetY = [(y / 35)**2 for y in np.linspace(99,30,70)]

spiralX = [(1-x/(2*np.pi))*np.cos(x+np.pi/4.) for x in np.linspace(0.,2*np.pi,180)] 
spiralY = [(1-y/(2*np.pi))*np.sin(y+np.pi/4.) for y in np.linspace(0.,2*np.pi,180)] 

for i in range(len(parabolaX)):
    fullPathX.append(parabolaX[i])
    fullPathY.append(parabolaY[i])
for i in range(len(parabolaRetX)):
    fullPathX.append(parabolaRetX[i])
    fullPathY.append(parabolaRetY[i])
for i in range(len(unitCircleX)):
    fullPathX.append(unitCircleX[i])
    fullPathY.append(unitCircleY[i])
for i in range(len(spiralX)):
    fullPathX.append(spiralX[i])
    fullPathY.append(spiralY[i])
'''
 ## Path 2 
spiralX = [2*(1-x/(2*np.pi))*np.cos(4*x) for x in np.linspace(0.,2*np.pi,360)] 
spiralY = [2*(1-y/(2*np.pi))*np.sin(4*y) for y in np.linspace(0.,2*np.pi,360)]
straightLineX = [2*(x/30) for x in np.linspace(0,29,30)]

for i in range(len(spiralX)):
    fullPathX.append(spiralX[i])
    fullPathY.append(spiralY[i])
for i in range(len(straightLineX)):
    fullPathX.append(straightLineX[i])
    fullPathY.append(0.)

### Path 3 

#for i in np.linspace(0,359,3600):
#    x,y = LJLocation(2.,3.,3*np.pi/4.,i) 
#    fullPathX.append(x)
#    fullPathY.append(y)


###
fig, ax = plt.subplots(figsize = (4,4))

ax.set_aspect(aspect='equal')
ax.set_xlabel('Real')
ax.set_ylabel('Imag')
ax.set_xlim(-2,2)
ax.set_ylim(-2,2)
ax.grid("on", lw = 0.5)
ax.axvline(x=0, ymin = - 1, ymax = 10, c ='black', lw = 1, ls = ':')
ax.axhline(y=0, xmin = - 1, xmax = 10, c ='black', lw = 1, ls = ':')


#ax.scatter(parabolaX,parabolaY, s= 0.1, c= 'r')
#ax.scatter(unitCircleX,unitCircleY, s = 0.1, c = 'purple')
#ax.scatter(fullPathX,fullPathY,s = 0.15, c = 'purple')


#seedLocation, = ax.plot(fullPathX[0],fullPathY[0], label= "c = ("+str(fullPathX[0])+","+str(fullPathY[0])+")", marker="o", ms=10, markerfacecolor="None",
#         markeredgecolor='red', markeredgewidth=1)

#ani = animation.FuncAnimation(fig, animate, repeat=True, frames=len(fullPathX)-1, interval = 20 )
#plt.tight_layout()
#ani.save(r"paths/SpiralInwards4Cycle2Radius360ImagesPath.gif",writer=writergif)
#plt.show()
#plt.clf()
embedGifs("SpiralInwards4Cycle2Radius360Images")