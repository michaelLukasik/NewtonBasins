import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import imageio
from PIL import Image, ImageSequence

writergif = animation.PillowWriter(fps=50)

def LJLocation(a,b,d,index, ph):
    return a*np.cos( 2*np.pi * (index/180) + d) , b*np.sin( 2*np.pi * (index/180))

def animate(i):
    x,y = pathLocation(i)
    seedLocation.set_data([x[i],y[i]])
    return seedLocation

def pathLocation (i):
    return np.array([fullPathX,fullPathY])

def examineGifs():
    gif1 = imageio.get_reader('Path1.gif')
    gif2 = imageio.get_reader(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\Path1GIF.gif")
    print(gif1.get_length())
    print(gif2.get_length())
    
    
def embedGifs():
    canvass = Image.open(r"BLANK.png")
    path   =  Image.open('Path1.gif')
    pattern =  Image.open(r"C:\Users\Michael\Documents\Programming\NewtonBasinImages\Bessel\variableOffset\Path1GIF.gif")
    totalH = max(path.size[0] ,pattern.size[0])
    totalW = path.size[1] + pattern.size[1]
    print(path.size[0],  pattern.size[0])
    print(path.size[1],  pattern.size[1])
    frames = []
    for i in range(530):
        path.seek(i % path.n_frames)
        pattern.seek(i % pattern.n_frames)
        frame = Image.new('RGB', (totalW-200, totalH), (255, 255, 255))
        frame.paste(path, (0,0), mask=path.convert("RGBA"))
        frame.paste(pattern, (path.size[1]-200,0), mask=pattern.convert("RGBA"))
        frames.append(frame)
    frames[0].save("testint.gif", save_all=True, duration=20, loop=0, append_images=frames[1:])


    
    


unitCircleX = [np.cos(x) for x in np.linspace(np.pi/4.,2*np.pi+np.pi/4.,180)] 
unitCircleY = [np.sin(y) for y in np.linspace(np.pi/4.,2*np.pi+np.pi/4.,180)] 

parabolaX = [(x / 35) for x in np.linspace(0,99,100)]
parabolaY = [(y / 35)**2 for y in np.linspace(0,99,100)]

parabolaRetX = [(x / 35) for x in np.linspace(99,30,70)]
parabolaRetY = [(y / 35)**2 for y in np.linspace(99,30,70)]

spiralX = [(1-x/(2*np.pi))*np.cos(x+np.pi/4.) for x in np.linspace(0.,2*np.pi,180)] 
spiralY = [(1-y/(2*np.pi))*np.sin(y+np.pi/4.) for y in np.linspace(0.,2*np.pi,180)] 

fullPathX = []
fullPathY = []

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



fig, ax = plt.subplots(figsize = (3,6))

ax.set_aspect(aspect='equal')
ax.set_xlabel('Real')
ax.set_ylabel('Imag')
ax.set_xlim(-1.1,3.5)
ax.set_ylim(-1.1,9.0)
ax.grid("on", lw = 0.5)

ax.axvline(x=0, ymin = - 1, ymax = 10, c ='black', lw = 1, ls = ':')
ax.axhline(y=0, xmin = - 1, xmax = 10, c ='black', lw = 1, ls = ':')


ax.scatter(parabolaX,parabolaY, s= 0.1, c= 'r')
ax.scatter(unitCircleX,unitCircleY, s = 0.1, c = 'purple')
ax.scatter(spiralX,spiralY,s = 0.1, c = 'pink')


#seedLocation, = ax.plot(fullPathX[0],fullPathY[0], marker="o", ms=10, markerfacecolor="None",
#         markeredgecolor='red', markeredgewidth=2)
#print(len(fullPathX)-1)
#ani = animation.FuncAnimation(fig, animate, repeat=True, frames=len(fullPathX)-1, interval = 20 )
#plt.tight_layout()
#ani.save("Path1.gif",writer=writergif)
#plt.show()
#plt.clf()
examineGifs()
embedGifs()