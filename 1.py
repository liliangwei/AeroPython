import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 50
xStart,xEnd = -2.0,2.0
yStart,yEnd = -1.0,1.0
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
print 'x= ',x
print 'y= ',y
X,Y = np.meshgrid(x,y)

#plotting the mesh

size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.scatter(X,Y,s=10,c='#CD2305',marker='o',linewidth=0.1)


strengthSource = 5.0
xSource, ySource =-1.0,0.0

uSource = np.empty((N,N),dtype=float)
vSource = np.empty((N,N),dtype=float)

#computing the velocity components at every point on the mesh grid
for i in range(N):
    for j in range(N):
        uSource[i,j] = strengthSource/(2*pi)\
            *(X[i,j]-xSource)/((X[i,j]-xSource)**2+(Y[i,j]-ySource)**2)
            
        vSource[i,j] = strengthSource/(2*pi)\
            *(Y[i,j]-ySource)/((X[i,j]-xSource)**2+(Y[i,j]-ySource)**2)

# plotting the streamlines
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uSource,vSource,density=2.0,linewidth=1,arrowsize=2,arrowstyle='->')
plt.scatter(xSource,ySource,c='#CD2305',s=80,marker='o')


strengthSink = -5.0
xSink,ySink = 1.0,0.0

uSink = np.empty((N,N),dtype=float)
vSink = np.empty((N,N),dtype=float)

for i in range(N):
    for j in range(N):
        uSink[i,j] = strengthSink/(2*pi)\
            *(X[i,j]-xSink)/((X[i,j]-xSink)**2+(Y[i,j]-ySink)**2)
        
        vSink[i,j] = strengthSink/(2*pi)\
            *(Y[i,j]-ySink)/((X[i,j]-xSink)**2+(Y[i,j]-ySink)**2)
            
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uSink,vSink,density=2.0,linewidth=1,arrowsize=2,arrowstyle='->')
plt.scatter(xSink,ySink,c='#CD2305',s=80,marker='o')

uPair = np.empty_like(uSource)
vPair = np.empty_like(vSource)

uPair = uSource+uSink
vPair = vSource+vSink

size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uPair,vPair,density=2.0,linewidth=1,arrowsize=2,arrowstyle='->')
plt.scatter([xSource,xSink],[ySource,ySink],c='#CD2305',s=80,marker='o')

 # challenge question
 # Contour plot of potential
   
potentialsource = np.empty_like(uPair)
potentialsink = np.empty_like(uPair)
potential = np.empty_like(uPair)

potentialsource = strengthSource/(2*pi)*np.log(np.sqrt((X-xSource)**2+(Y-ySource)**2))
potentialsink = strengthSink/(2*pi)*np.log(np.sqrt((X-xSink)**2+(Y-ySink)**2))
potential = potentialsource + potentialsink


plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.contourf(X,Y,potential,levels=np.linspace(-2.0,2.0,100),extend='both')
plt.scatter([xSource,xSink],[ySource,ySink],c='black',s=80,marker='o')
            
plt.show()