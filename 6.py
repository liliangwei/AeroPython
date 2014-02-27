import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 50
xStart,xEnd = -2.0,2.0
yStart,yEnd = -1.0,1.0
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)

kappa =1.0
xDoublet,yDoublet = 0.0,0.0

def getVelocityDoublet(strength,xd,yd,X,Y):
    u = -strength/(2*pi)*((X-xd)**2-(Y-yd)**2)/((X-xd)**2+(Y-yd)**2)**2
    v = -strength/(2*pi)*2*(X-xd)*(Y-yd)/((X-xd)**2+(Y-yd)**2)**2
    return u,v

def getStreamFunctionDoublet(strength,xd,yd,X,Y):
    psi = -strength/(2*pi)*(Y-yd)/((X-xd)**2+(Y-yd)**2)
    return psi

uDoublet,vDoublet = getVelocityDoublet(kappa,xDoublet,yDoublet,X,Y)

psiDoublet = getStreamFunctionDoublet(kappa,xDoublet,yDoublet,X,Y)

Uinf = 1.0

uFreestream = Uinf*np.ones((N,N),dtype=float)
vFreestream = np.zeros((N,N),dtype = float)

psiFreestream = Uinf*Y

u = uFreestream+uDoublet
v = vFreestream+vDoublet
psi = psiFreestream+psiDoublet

size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xDoublet,yDoublet,c='yellow',s=80,marker='o')

#cylinder radius
R = sqrt(kappa/(2*pi*Uinf))
circle = plt.Circle((0,0),radius=R,color='r',alpha=0.5)
plt.gca().add_patch(circle)

#stagnation points
xStagn1,yStagn1 = sqrt(kappa/(2*pi*Uinf)),0
xStagn2,yStagn2 = -sqrt(kappa/(2*pi*Uinf)),0
plt.scatter([xStagn1,xStagn2],[yStagn1,yStagn2],c='black',s=80,marker='o');

#Vortex
gamma = 3.0
xVortex,yVortex = 0.0,0.0

def getVelocityVortex(strength,xv,yv,X,Y):
    u = strength/(2*pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = -strength/(2*pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    return u,v
    
def getStreamFunctionVortex(strength,xv,yv,X,Y):
    psi = strength/(4*pi)*np.log((X-xv)**2+(Y-yv)**2)
    return psi
    
uVortex,vVortex = getVelocityVortex(gamma,xVortex,yVortex,X,Y)

psiVortex = getStreamFunctionVortex(gamma,xVortex,yVortex,X,Y)

u = uFreestream+uDoublet+uVortex
v = vFreestream+vDoublet+vVortex
psi = psiFreestream+psiDoublet+psiVortex

#cyliner radius
R = sqrt(kappa/(2*pi*Uinf))

#stagnation points
xStagn1,yStagn1 = sqrt(R**2-(gamma/(4*pi*Uinf))**2),-gamma/(4*pi*Uinf)
xStagn2,yStagn2 = -sqrt(R**2-(gamma/(4*pi*Uinf))**2),-gamma/(4*pi*Uinf)

size=10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
circle = plt.Circle((0,0),radius=R,color='r',alpha=0.5)
plt.gca().add_patch(circle)
plt.scatter(xVortex,yVortex,c='yellow',s=80,marker='o')
plt.scatter([xStagn1,xStagn2],[yStagn1,yStagn2],c='black',s=80,marker='o')

#Pressure coefficiet
theta = np.linspace(0,2*pi,100)
utheta = -2*Uinf*np.sin(theta)-gamma/(2*pi*R)

Cp = 1-(utheta/Uinf)**2

# if there was no vortex
utheta_noVortex = -2*Uinf*np.sin(theta)
Cp_noVortex = 1-(utheta_noVortex/Uinf)**2

#plotting
size=6
plt.figure(figsize=(size,size))
plt.grid(True)
plt.xlabel(r'$\theta$',fontsize=18)
plt.ylabel(r'$C_p$',fontsize=18)
plt.xlim(theta.min(),theta.max())
plt.plot(theta,Cp,color='r',linewidth=2,linestyle='-')
plt.plot(theta,Cp_noVortex,color='blue',linewidth=2,linestyle='--')
plt.legend(['with vortex','without vortex'],loc='best',prop={'size':16})

# Lift and Drag

#Without vertex
#drag = -(1./3.)*np.sin(3.*theta)*0.5*R*Uinf**2
#lift = (1./3.)*(np.cos(3.*theta)-6.*np.cos(theta))*0.5*R*Uinf**2

#With vortex
x_c = gamma/2/pi/R

drag = 4/3*(np.sin(theta))**3+0.2*np.sin(theta)-x_c*(np.cos(theta))**2*R*.5
lift = (x_c*theta-.5*np.sin(2.*theta)-3.2*np.cos(theta)+1/3*np.cos(3.*theta))*R*.5

size=6
plt.figure(figsize=(size,size))
plt.grid(True)
plt.xlabel(r'$\theta$',fontsize=18)
plt.ylabel(r'$Force$',fontsize=18)
plt.xlim(theta.min(),theta.max())
plt.plot(theta,drag,color='r',linewidth=2,linestyle='-')
plt.plot(theta,lift,color='blue',linewidth=2,linestyle='--')
plt.legend(['Drag','lift'],loc='best',prop={'size':16})

plt.show()
