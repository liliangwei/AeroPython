# Vortex

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
xVortex,yVortex = 0.0,0.0

def getVelocityVortex(strength,xv,yv,X,Y):
    u = strength/(2*pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v =-strength/(2*pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    return u,v

def getStreamFunctionVortex(strength,xv,yv,X,Y):
    psi = strength/(4*pi)*np.log((X-xv)**2+(Y-yv)**2)
    return psi

uVortex,vVortex = getVelocityVortex(gamma,xVortex,yVortex,X,Y)

psiVortex = getStreamFunctionVortex(gamma,xVortex,yVortex,X,Y)

#plotting
size=10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uVortex,vVortex,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xVortex,yVortex,c='yellow',s=80,marker='o');

strengthSink = -1.0
xSink,ySink = 0.0,0.0

def getVelocitySink(strength,xs,ys,X,Y):
    u = strength/(2*pi)*(X-xs)/((X-xs)**2+(Y-ys)**2)
    v = strength/(2*pi)*(Y-ys)/((X-xs)**2+(Y-ys)**2)
    return u,v

def getStreamFunctionSink(strength,xs,ys,X,Y):
    psi = strength/(2*pi)*np.arctan2((Y-ys),(X-xs))
    return psi

uSink,vSink = getVelocitySink(strengthSink,xSink,ySink,X,Y)

psiSink = getStreamFunctionSink(strengthSink,xSink,ySink,X,Y)


u = uVortex+uSink
v = vVortex+vSink

psi = psiVortex+psiSink

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xVortex,yVortex,c='black',s=80,marker='o');

#plt.show()


