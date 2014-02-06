import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 200
xStart,xEnd = -4.0,4.0
yStart,yEnd = -2.0,2.0
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)

np.shape(X)

Uinf = 1.0
alpha = 0.0

uFreestream = Uinf*cos(alpha*pi/180)*np.ones((N,N),dtype=float)
vFreestream = Uinf*sin(alpha*pi/180)*np.ones((N,N),dtype=float)

psiFreestream = + Uinf*cos(alpha*pi/180)*Y\
                - Uinf*sin(alpha*pi/180)*X
                
# function to compute the velocity field of a source/sink
def getVelocity(strength,xs,ys,X,Y):
    u = strength/(2*pi)*(X-xs)/((X-xs)**2+(Y-ys)**2)
    v = strength/(2*pi)*(Y-ys)/((X-xs)**2+(Y-ys)**2)
    return u,v
    
def getStreamFunction(strength,xs,ys,X,Y):
    psi = strength/(2*pi)*np.arctan2((Y-ys),(X-xs))
    return psi

strengthSource = 5.0
xSource,ySource = -1.0,0.0

uSource,vSource = getVelocity(strengthSource,xSource,ySource,X,Y)

psiSource = getStreamFunction(strengthSource,xSource,ySource,X,Y)

# superposition of the source on the freestream

u = uFreestream + uSource
v = vFreestream + vSource
psi = psiFreestream + psiSource

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xSource,ySource,c='#CD2305',s=80,marker='o')

# computing the stagnation point
xStagnation = xSource - strengthSource/(2*pi*Uinf)*cos(alpha*pi/180)
yStagnation = ySource - strengthSource/(2*pi*Uinf)*sin(alpha*pi/180)

# adding the stagnation point to the figure
plt.scatter(xStagnation,yStagnation,c='yellow',s=80,marker='o')

# adding the dividing line to the figure
if (alpha==0.0):
    plt.contour(X,Y,psi,\
            levels=[-strengthSource/2,+strengthSource/2],\
            colors='#CD2305',linewidths=2,linestyles='--')
            
# Sink Part
strengthSink = -5.0
xSink,ySink = 1.0,0.0

uSink,vSink = getVelocity(strengthSink,xSink,ySink,X,Y)

psiSink = getStreamFunction(strengthSink,xSink,ySink,X,Y)


u = uFreestream + uSource + uSink
v = vFreestream + vSource + vSink

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel

plt.show()