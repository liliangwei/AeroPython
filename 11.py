# Source-vortex panel method

import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt

# read of the geometry
coords = np.loadtxt(fname='C:/Users/llwei89/Documents/Github/AeroPython/resources/n0012.dat')
xp,yp = coords[:,0],coords[:,1]

# plotting the geometry
valX,valY = 0.1,0.2
xmin,xmax = min(xp),max(xp)
ymin,ymax = min(yp),max(yp)
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)

size=10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2);

# class Panel containing the info about one panel
class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya
        self.xb,self.yb = xb,yb
        self.xc,self.yc = (xa+xb)/2.,(ya+yb)/2.
        self.length = sqrt((xb-xa)**2+(yb-ya)**2)
        
        # orientation of panel
        if (xb-xa<=0.) : self.beta = acos((yb-ya)/self.length)
        elif (xb-xa>0.) : self.beta = pi+acos(-(yb-ya)/self

plt.show()