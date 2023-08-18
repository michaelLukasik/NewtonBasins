import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import comb

def bernstein_poly(i, n, t):
    return comb(n, i) * ( t**(n-i) ) * (1 - t)**i


def bezier_curve(points, nTimes):

    nPoints = len(points)
    xPoints = np.array([p[0] for p in points])
    yPoints = np.array([p[1] for p in points])

    t = np.linspace(0.0, 1.0, nTimes)

    polynomial_array = np.array([ bernstein_poly(i, nPoints-1, t) for i in range(0, nPoints)   ])

    xvals = np.dot(xPoints, polynomial_array)
    yvals = np.dot(yPoints, polynomial_array)
    return np.flip(xvals), np.flip(yvals)

class Barcalona:
    def __init__(self, offsetStartPoint, rotation):
        self.length = 	4.657 # km
        self.turns = 14 ## No chicane at turn 14
        self.nBezierPoints = 10
        self.FinishLine = np.array([[75.,0.]]) + offsetStartPoint * rotation
        
        
class sector(): 
        def __init__(self):
            self.start = np.array([[0.,0.]])
            self.end = np.array([[0.,0.]])
            
            self.sectorValuesX = np.array([[]])
            self.sectorValuesY = np.array([[]])

        
        def getSectorValues(self):
            return self.sectorValuesX, self.sectorValuesY
            
class section(): 
        
    def __init__(self):
        self.start = np.array([[0.,0.]])
        self.end = np.array([[1.,0.]])
        self.controlPoints = np.array([[0.5,0.]])
        self.nPointsAlongSection = 20
        self.sectionValuesX = np.array([[]])
        self.sectionValuesY = np.array([[]])

    
    def getBezierPoints(self):
        return np.concatenate((self.start,self.controlPoints,self.end), axis =0)
    
    def getSectionValues(self, sectorObject):
        points = self.getBezierPoints()
        self.sectionValuesX , self.sectionValuesY  = bezier_curve(points, self.nPointsAlongSection)
        sectorObject.sectorValuesX = np.append(sectorObject.sectorValuesX, self.sectionValuesX[:-1])
        sectorObject.sectorValuesY = np.append(sectorObject.sectorValuesY, self.sectionValuesY[:-1])
        return self.sectionValuesX , self.sectionValuesY
        
            
        
class turn(): 
    def __init__(self):
        self.start = np.array([[0.,0.]])
        self.end = np.array([[1.,0.]])
        self.controlPoints = np.array([[0.5,0.]])
        self.nPointsAlongTurn = 20
        self.turnValuesX = np.array([[]])
        self.turnValuesY = np.array([[]])
        
        
    
    def getBezierPoints(self):
        return np.concatenate((self.start,self.controlPoints,self.end), axis =0)
    
    def getTurnValues(self, sectorObject):
        points = self.getBezierPoints()
        self.turnValuesX , self.turnValuesY  = bezier_curve(points, self.nPointsAlongTurn)
        sectorObject.sectorValuesX = np.append(sectorObject.sectorValuesX, self.turnValuesX[:-1], )
        sectorObject.sectorValuesY = np.append(sectorObject.sectorValuesY, self.turnValuesY[:-1], )
        return self.turnValuesX , self.turnValuesY
    
    



def plotBySector():
    Sector1X,Sector1Y = Sector1.getSectorValues()
    Sector2X,Sector2Y = Sector2.getSectorValues()
    Sector3X,Sector3Y = Sector3.getSectorValues()
    fig, ax = plt.subplots ()
    ax.scatter(Sector1X[0],Sector1Y[0], s = 500, c = 'purple', marker = "|")
    ax.scatter(Sector1X,Sector1Y, s = 3, c = 'r')
    ax.scatter(Sector2X,Sector2Y, s = 3, c = 'b')
    ax.scatter(Sector3X,Sector3Y, s = 3, c = 'y')
    ax.set_xlim(0,122)
    ax.set_ylim(-3,37)
    ax.grid(1)
    ax.set_aspect(1)
    plt.title("Circuit de Barcelona-Catalunya")
    plt.tight_layout()
    plt.savefig(r"CircuitdeBarcelona-Catalunya.png")
    plt.show()
    
def exportTrackPoints():
    x = np.concatenate([Sector1X, Sector2X, Sector3X], axis = 0)
    y = np.concatenate([Sector1Y, Sector2Y, Sector3Y], axis = 0)
    df = pd.DataFrame(y,x)
    df.to_csv("barcalonaTrackPoints.csv", header = None)    
                
#################################################################


offset = np.array([[0.,0.]])
rotation = np.array([[0.,0.]])  
BarcalonaTrack = Barcalona(offset,rotation)

##########################    SECTOR 1  ########################## We drive CLOCKWISE
Sector1 = sector()

### Section 0

s0 = section()
s0.start =  BarcalonaTrack.FinishLine
s0.controlPoints =  np.array([[55.,0.]])
s0.end =  np.array([[25.,0.0]])
s0.nPointsAlongSection = 20
s0.getBezierPoints()
s0.getSectionValues(Sector1)

print("1     ,   ", Sector1.sectorValuesX)
### Turn 1 

t1 = turn()
t1.start = s0.end
t1.end = np.array([[20.,5.]])
t1.controlPoints =  np.array([[20.,0.]])
t1.nPointsAlongTurn = 20
t1.getBezierPoints()
t1.getTurnValues(Sector1)

## Turn 2 

t2 = turn()
t2.start = t1.end
t2.end = np.array([[12.,12.]])
t2.controlPoints =  np.array([[20.,10.]])
t2.nPointsAlongTurn = 10
t2.getBezierPoints()
t2.getTurnValues(Sector1)


## Turn 3

t3 = turn()
t3.start = t2.end
t3.end = np.array([[12.,32.]])
t3.controlPoints =  np.array([[2.,25.]])
t3.nPointsAlongTurn = 15
t3.getBezierPoints()
t3.getTurnValues(Sector1)


## Section 3a

s3a = section()
s3a.start = t3.end
s3a.end = np.array([[35.,32.]])
s3a.controlPoints =  np.array([[24.,34.5]])
s3a.nPointsAlongSection = 15
s3a.getBezierPoints()
s3a.getSectionValues(Sector1)
 
 
##########################    SECTOR 2  ########################## 
Sector2 = sector()

## Turn 4

t4 = turn()
t4.start = s3a.end
t4.end = np.array([[33.5,20.]])
t4.controlPoints =  np.array([[45.,28.5]])
t4.nPointsAlongTurn = 10
t4.getBezierPoints()
t4.getTurnValues(Sector2)

## Section 4a

s4a = section()
s4a.start = t4.end
s4a.end = np.array([[20.5,19.]])
s4a.controlPoints =  np.array([[25.,19.5]])
s4a.nPointsAlongSection = 10
s4a.getBezierPoints()
s4a.getSectionValues(Sector2)



## Turn 5

t5 = turn()
t5.start = s4a.end
t5.end = np.array([[20.0,16.]])
t5.controlPoints =  np.array([[15.,18.5]])
t5.nPointsAlongTurn = 10
t5.getBezierPoints()
t5.getTurnValues(Sector2)

## Section 5a

s5a = section()
s5a.start = t5.end
s5a.end = np.array([[30.,10.]])
s5a.controlPoints =  np.array([[25.,13]])
s5a.nPointsAlongSection = 10
s5a.getBezierPoints()
s5a.getSectionValues(Sector2)


## Turn 6
t6 = turn()
t6.start = s5a.end
t6.end = np.array([[45.0,9.0]])
t6.controlPoints =  np.array([[32.,8.5]])
t6.nPointsAlongTurn = 10
t6.getBezierPoints()
t6.getTurnValues(Sector2)

## Turn 7 

t7 = turn()
t7.start = t6.end
t7.end = np.array([[48.0,14.0]])
t7.controlPoints =  np.array([[49.,11.5]])
t7.nPointsAlongTurn = 7
t7.getBezierPoints()
t7.getTurnValues(Sector2)


## Turn 8

t8 = turn()
t8.start = t7.end
t8.end = np.array([[55.0,27.0]])
t8.controlPoints =  np.array([[47.,15.5]])
t8.nPointsAlongTurn = 10
t8.getBezierPoints()
t8.getTurnValues(Sector2)


## Turn 9

t9 = turn()
t9.start = t8.end
t9.end = np.array([[72.0,23.0]])
t9.controlPoints =  np.array([[62.5,35.5]])
t9.nPointsAlongTurn = 10
t9.getBezierPoints()
t9.getTurnValues(Sector2)


## Section 9a

s9a = turn()
s9a.start = t9.end
s9a.end = np.array([[95.0,11.0]])
s9a.controlPoints =  np.array([[90.0,12.2]])
s9a.nPointsAlongTurn = 12
s9a.getBezierPoints()
s9a.getTurnValues(Sector2)


##########################    SECTOR 3  ########################## 
Sector3 = sector()

## Turn 10

t10 = turn()
t10.start = s9a.end
t10.end = np.array([[104.0,16.0]])
t10.controlPoints =  np.array([[104.0,8.5]])
t10.nPointsAlongTurn = 7
t10.getBezierPoints()
t10.getTurnValues(Sector3)

## Turn 11

t11 = turn()
t11.start = t10.end
t11.end = np.array([[98.0,21.5]])
t11.controlPoints =  np.array([[103.0,21.]])
t11.nPointsAlongTurn = 7
t11.getBezierPoints()
t11.getTurnValues(Sector3)


## Turn 12

t12 = turn()
t12.start = t11.end
t12.end = np.array([[92.5,30.]])
t12.controlPoints =  np.array([[90,23.]])
t12.nPointsAlongTurn = 12
t12.getBezierPoints()
t12.getTurnValues(Sector3)


## section 12a

s12a = turn()
s12a.start = t12.end
s12a.end = np.array([[115.0,30.]])
s12a.controlPoints =  np.array([[97.5,40.5]])
s12a.nPointsAlongTurn = 15
s12a.getBezierPoints()
s12a.getTurnValues(Sector3)


## Turn 13

t13 = turn()
t13.start = s12a.end
t13.end = np.array([[120.,12.]])
t13.controlPoints =  np.array([[121,27.]])
t13.nPointsAlongTurn = 12
t13.getBezierPoints()
t13.getTurnValues(Sector3)


## Turn 14

t14 = turn()
t14.start = t13.end
t14.end =  np.array([[95.,0.]])
t14.controlPoints =  np.array([[122,-1.]])
t14.nPointsAlongTurn = 20
t14.getBezierPoints()
t14.getTurnValues(Sector3)


## section 14a

s14a = turn()
s14a.start = t14.end
s14a.end =  BarcalonaTrack.FinishLine
s14a.controlPoints =  np.array([[90.,0.]])
s14a.nPointsAlongTurn = 8
s14a.getBezierPoints()
s14a.getTurnValues(Sector3)



Sector1X,Sector1Y = Sector1.getSectorValues()
Sector2X,Sector2Y = Sector2.getSectorValues()
Sector3X,Sector3Y = Sector3.getSectorValues()
plotBySector()
#exportTrackPoints()



