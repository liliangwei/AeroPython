import numpy as np
import matplotlib.pyplot as plt
from math import *
from IPython.core.display import clear_output

N = 50
xStart,xEnd = -2.0,2.0
yStart,yEnd = -1.0,1.0
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)

# Class of Source
class Source:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
#get velocity field
    def velocity(self,x,y):
        self.u = self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
        self.v = self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
#get streamfunction
    def streamFunction(self,x,y):
        self.psi = self.strength/(2*pi)*np.arctan2((Y-self.y),(X-self.x))

strengthSource = 1.0
xSource,ySource = 0.0,0.5
source = Source(strengthSource,xSource,ySource)
source.velocity(X,Y)
source.streamFunction(X,Y)

#Source Image
sourceImage = Source(strengthSource,xSource,-ySource)
sourceImage.velocity(X,Y)
sourceImage.streamFunction(X,Y)

#Superposition of velocity
u = source.u + sourceImage.u
v = source.v + sourceImage.v
psi = source.psi + sourceImage.psi

# First plotting

size = 10
plt.figure(num=0,figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(source.x,source.y,c='r',s=80,marker='o')
plt.scatter(sourceImage.x,sourceImage.y,c='r',s=80,marker='D')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4);

#Class of Vortex
class Vortex:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
    #get velocity field
    def velocity(self,X,Y):
        self.u = self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
        self.v = -self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
    
    #get streamline
    def streamFunction(self,X,Y):
        self.psi = -self.strength/(4*pi)*np.log((X-self.x)**2+(Y-self.y)**2)
    
strengthVortex = 1.0
xVortex,yVortex=0.0,0.5
vortex = Vortex(strengthVortex,xVortex,yVortex)
vortex.velocity(X,Y)
vortex.streamFunction(X,Y)

vortexImage = Vortex(-strengthVortex,xVortex,-yVortex)
vortexImage.velocity(X,Y)
vortexImage.streamFunction(X,Y)

#superposition of vortex
u = vortex.u + vortexImage.u
v = vortex.v + vortexImage.v
psi = vortex.psi + vortexImage.psi

#plotting
size = 10
plt.figure(num=1,figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(vortex.x,vortex.y,c='r',s=80,marker='o')
plt.scatter(vortexImage.x,vortexImage.y,c='r',s=80,marker='D')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4);

#Vortex pair
xVortex1,xVortex2 = -0.1,0.1
yVortex1,yVortex2 = 0.5,0.5

vortex1 = Vortex(strengthVortex,xVortex1,yVortex1)
vortex2 = Vortex(-strengthVortex,xVortex2,yVortex2)
vortex1.velocity(X,Y)
vortex1.streamFunction(X,Y)
vortex2.velocity(X,Y)
vortex2.streamFunction(X,Y)

vortexImage1 = Vortex(-strengthVortex,xVortex1,-yVortex1)
vortexImage2 = Vortex(strengthVortex,xVortex2,-yVortex2)
vortexImage1.velocity(X,Y)
vortexImage1.streamFunction(X,Y)
vortexImage2.velocity(X,Y)
vortexImage2.streamFunction(X,Y)

#superposition pair vortex
u = vortex1.u + vortex2.u + vortexImage1.u + vortexImage2.u
v = vortex1.v + vortex2.v + vortexImage1.v + vortexImage2.v
psi = vortex1.psi + vortex2.psi + vortexImage1.psi + vortexImage2.psi

#plotting
size = 10
plt.figure(num=2,figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(vortex1.x,vortex1.y,c='r',s=80,marker='o')
plt.scatter(vortex2.x,vortex2.y,c='g',s=80,marker='o')
plt.scatter(vortexImage1.x,vortexImage1.y,c='r',s=80,marker='D')
plt.scatter(vortexImage2.x,vortexImage2.y,c='g',s=80,marker='D')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4)
plt.axvline(0.0,color='k',linestyle='--',linewidth=2);

#Doublet near a plane
Uinf = 1.0
uFreestream = Uinf*np.ones((N,N),dtype=float)
vFreestream = np.zeros((N,N),dtype=float)
psiFreestream = Uinf*Y

# Doublet class
class Doublet:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
        
    #Velocity Field
    def velocity(self,X,Y):
        self.u = -self.strength/(2*pi)*((X-self.x)**2-(Y-self.y)**2)/((X-self.x)**2+(Y-self.y)**2)
        self.v = -self.strength/(2*pi)*2*(X-self.x)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)**2
        
    def streamFunction(self,X,Y):
        self.psi = -self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
    

strengthDoublet = 1.0
xDoublet,yDoublet = 0.0,0.3
doublet = Doublet(strengthDoublet,xDoublet,yDoublet)
doublet.velocity(X,Y)
doublet.streamFunction(X,Y)

doubletImage = Doublet(strengthDoublet,xDoublet,-yDoublet)
doubletImage.velocity(X,Y)
doubletImage.streamFunction(X,Y)

#Superposition  of the doublet
u = uFreestream + doublet.u + doubletImage.u
v = vFreestream + doublet.v + doubletImage.v
psi = psiFreestream + doublet.psi + doubletImage.psi

#plotting
size = 10
plt.figure(num=3,figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(doublet.x,doublet.y,c='r',s=80,marker='o')
plt.scatter(doubletImage.x,doubletImage.y,c='r',s=80,marker='x')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4);

        

plt.show()