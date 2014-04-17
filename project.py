# Project 
# Starson
# Two element panel method

import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt

# read of the geometry
coords = np.loadtxt(fname='C:/Users/llwei89/Documents/Github/AeroPython/resources/n0012.dat')
#coords = np.loadtxt(fname='/home/starson/AeroPython/resources/n0012.dat')
xp1,yp1 = coords[:,0],coords[:,1]            # read in the original airfoil

# Scale the second flap
xp2 = 0.5*np.copy(xp1)                       # scaled x coordinate of flap
yp2 = 0.5*np.copy(yp1)                       # scaled y coordinate of flap

# rotate the flap angle delta 35 degree down
delta = 35

x1 = np.zeros_like(xp2)      
y1 = np.zeros_like(yp2)

L = len(xp1)

for i in range(L):
    r  = sqrt(xp2[i]**2+yp2[i]**2)
    thetatr = np.arctan2(yp2[i],xp2[i])
    theta = thetatr - delta*pi/180.0
    x1[i] = r*cos(theta)
    y1[i] = r*sin(theta)

xp2 = x1
yp2 = y1

xp2 = xp2 + 0.92                            # translated x coordinate
yp2 = yp2 -0.06                             # translated y coordinate

# plotting the geometry
valX,valY = 0.2,0.4                          #value for plot margin    
xmin,xmax = min(xp1),max(xp2)
ymin,ymax = min(yp2),max(yp1)
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)

size=10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp1,yp1,xp2,yp2,'k-',linewidth=2);



plt.show()