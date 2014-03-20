# Souece Sheet

import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy import integrate

N = 100
xStart,xEnd = -1.0,1.0
yStart,yEnd = -1.5,1.5

x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)

Uinf = 1.0

uFreestream = Uinf*np.ones((N,N),dtype=float)
vFreestream = np.zeros((N,N),dtype=float)

# Class of Source
class Source:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
    #get velocity
    def velocity(self,X,Y):
        self.u = self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
        self.v = self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
    
    #get the stream-function
    def streamFunction(self,X,Y):
        self.psi = self.strength/(2*pi)*np.arctan2((Y-self.y),(X-self.x))

# definition of my sources
NSource = 11
strength = 5.0
strengthSource = strength/NSource
xSource = 0.0
ySource = np.linspace(-1.0,1.0,NSource)

sourceSheet = np.empty(NSource,dtype=object)
for i in range(NSource):
    sourceSheet[i] = Source(strengthSource,xSource,ySource[i])
    sourceSheet[i].velocity(X,Y)

# superposition
u = uFreestream.copy()
v = vFreestream.copy()

for i in sourceSheet:
    u = np.add(u,i.u)
    v = np.add(v,i.v)

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=3,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xSource*np.ones(NSource,dtype=float),ySource,c='r',s=80,marker='D')
cont = plt.contourf(X,Y,np.sqrt(u**2+v**2),levels=np.linspace(0.0,0.1,10))
cbar = plt.colorbar(cont)
cbar.set_label('U',fontsize=16);

a = 1.0
f = lambda x : x+a
print 'x + a =',f(1)

print integrate.quad(lambda x : x,0.0,1.0)[0]

sigma = 2.5     # strength of the source sheet

uPanel = np.empty((N,N),dtype=float)
vPanel = np.empty((N,N),dtype=float)

# boundaries of the sheet

ymin,ymax = -1.0,1.0

# computing the velocity field
for i in range(N):
    for j in range(N):
        
        func = lambda s : X[i,j]/(X[i,j]**2+(Y[i,j]-s)**2)
        uPanel[i,j] = sigma/(2*pi)*integrate.quad(func,ymin,ymax)[0]
        
        func = lambda s:(Y[i,j]-s)/(X[i,j]**2+(Y[i,j]-s)**2)
        vPanel[i,j] = sigma/(2*pi)*integrate.quad(func,ymin,ymax)[0]
        
u2 = np.add(uFreestream,uPanel)
v2 = np.add(vFreestream,vPanel)

#plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=3,linewidth=1,arrowsize=1,arrowstyle='->')
plt.axvline(0.0,(ymin-yStart)/(yEnd-yStart),(ymax-yStart)/(yEnd-yStart),c='r',linewidth=4)
cont = plt.contourf(X,Y,np.sqrt(u2**2+v2**2),levels= np.linspace(0.0,2.0,10))
cbar = plt.colorbar(cont)
cbar.set_label('U',fontsize=16);

plt.show()