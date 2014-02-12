import numpy as np
import matplotlib.pyplot as plt
from math import *

#from IPython.display import Image
#Image(filename='../resources/doubleSketch2.png')

N = 50
xStart,xEnd = -2.0,2.0
yStart,yEnd = -1.0,1.0
x = np.linspace(xStart,xEnd,N)
y = np.linspace(xStart,yEnd,N)
X,Y = np.meshgrid(x,y)

kappa = 1.0
xDoublet,yDoublet = 0.0,0.0

def getVelocityDoublet(strength,xd,yd,X,Y):
    u = -strength/(2*pi)*((X-xd)**2-(Y-yd)**2)/((X-xd)**2+(Y-yd)**2)**2
    v = -strength/(2*pi)*2*(X-xd)*(Y-yd)/((X-xd)**2+(Y-yd)**2)**2
    return u,v

def getStreamFunctionDoublet(strength,xd,yd,X,Y):
    psi = - strength/(2*pi)*(Y-yd)/((X-xd)**2+(Y-yd)**2)
    return psi
    
uDoublet,vDoublet = getVelocityDoublet(kappa,xDoublet,yDoublet,X,Y)

psiDoublet = getStreamFunctionDoublet(kappa,xDoublet,yDoublet,X,Y)

#plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)

plt.streamplot(X,Y,uDoublet,vDoublet,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xDoublet,yDoublet,c='yellow',s=80,marker='o')

Uinf = 1.0

uFreestream = Uinf*np.ones((N,N),dtype=float)
vFreestream = Uinf*np.zeros((N,N),dtype=float)

psiFreestream = Uinf*Y

# superimposition of the doublet on the freestream flow
u = uFreestream + uDoublet
v = vFreestream + vDoublet

psi = psiDoublet + psiFreestream

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.contour(X,Y,psi,levels=[0.0],colors='#CD2305',linewidths=2,linestyles='--')
plt.scatter(xDoublet,yDoublet,c='yellow',s=80,marker='o')

#stagnation points
xStagn1,yStagn1 = sqrt(kappa/(2*pi*Uinf)),0
xStagn2,yStagn2 = -sqrt(kappa/(2*pi*Uinf)),0
plt.scatter([xStagn1,xStagn2],[yStagn1,yStagn2],c='black',s=80,marker='o');

#compute the pressure coefficient

Cp = 1.0-(u**2+v**2)/Uinf**2

#plotting
size = 10
plt.figure(num=0,figsize=(1.1*size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
contf = plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
cbar = plt.colorbar(contf)
cbar.set_label(r'$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
plt.scatter(xDoublet,yDoublet,c='yellow',s=80,marker='o')
plt.contour(X,Y,psi,levels=[0,0],color='#CD2305',linewidths=2,linestyle='--')
plt.scatter([xStagn1,xStagn2],[yStagn1,yStagn2],c='black',s=80,marker='o');

plt.show()
