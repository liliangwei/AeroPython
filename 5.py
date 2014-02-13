# infinite row of vortices
# Starson

import numpy as np
import matplotlib.pyplot as plt
from math import *



N = 50
xStart,xEnd = -2.0,2.0
yStart,yEnd = -1.0,1.0
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)

gamma = 5.0
vnumber = 11
length = (xEnd - xStart)/vnumber
#------------------------finite number of vortex-------------------------------


xvortex = np.linspace(xStart,xEnd,vnumber)
yvortex = np.zeros(vnumber)

def getVelocityVortex(strength,xv,yv,X,Y):
    u = strength/(2*pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v =-strength/(2*pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    return u,v

def getStreamFunctionVortex(strength,xv,yv,X,Y):
    psi = strength/(4*pi)*np.log((X-xv)**2+(Y-yv)**2)
    return psi

u = np.zeros((N,N),dtype=float)
v = np.zeros((N,N),dtype=float)
psi = np.zeros((N,N),dtype=float)
for i in range(vnumber):
    uVortex,vVortex = getVelocityVortex(gamma,xvortex[i],yvortex[i],X,Y)
    psiVortex = getStreamFunctionVortex(gamma,xvortex[i],yvortex[i],X,Y) 
    u += uVortex
    v += vVortex
    psi += psiVortex


# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
for i in range(vnumber):
    plt.scatter(xvortex[i],yvortex[i],c='yellow',s=80,marker='o');
    
#---------infinite number of vortex-------------------------------------------


def getInfiniteVelocity(strength,a,X,Y):
    u = strength/(2*pi)*np.sinh(2*pi*Y/a)/(np.cosh(2*pi*Y/a)-np.cos(2*pi*X/a))
    v = -strength/(2*pi)*np.sin(2*pi*Y/a)/(np.cosh(2*pi*Y/a)-np.cos(2*pi*X/a))
    return u,v


inuVortex,invVortex = getInfiniteVelocity(gamma,length,X,Y)
    
# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,inuVortex,invVortex,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
for i in range(vnumber):
    plt.scatter(xvortex[i],yvortex[i],c='yellow',s=80,marker='o');

    

plt.show()
    
