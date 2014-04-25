# Project 
# Starson
# Two element panel method

import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt



#Set up free stream flow
class Freestream:
    def __init__(self,Uinf,alpha):
        self.Uinf = Uinf
        self.alpha = alpha*pi/180.0

Uinf = 5.0
alpha = 0       # the angle of attack from main elemnt chord line
AL = alpha / 57.2958    # get into radians


freestream = Freestream(Uinf,alpha)
# read of the geometry
#coords = np.loadtxt(fname='C:/Users/llwei89/Documents/Github/AeroPython/resources/n0012.dat')
coords = np.loadtxt(fname='/home/starson/AeroPython/resources/n0012.dat')
xp,yp = coords[:,0],coords[:,1]            # read in the original airfoil

# plotting the geometry
valX,valY = 0.2,0.4                          #value for plot margin    
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
        self.TH = atan2((yb-ya),(xb-xa))
        
        # orientation of panel
        if (xb-xa<=0.) : self.beta = acos((yb-ya)/self.length)
        elif (xb-xa>0.) : self.beta = pi+acos(-(yb-ya)/self.length)
        
        #location of the panel
        if (self.beta<=pi): self.loc = 'extrados'
        else: self.loc = 'intrados'
        
        self.gamma = 0.0                            # vortex strength
        self.vt = 0.0                               # tengential velocity
        self.Cp = 0.0                               # pressure coefficient
        
def definePanels(N,xp,yp):
    
    R = (max(xp)-min(xp))/2.0                       # radius of circle
    xCenter = (max(xp)+min(xp))/2.0                 # x center
    xCircle = xCenter + R*np.cos(np.linspace(0,2*pi,N+1))
    
    x = np.copy(xCircle)      # projection of the x-coord on the surface
    y = np.empty_like(x)      # initialization of the y-coord Numpy array
    
    xp,yp = np.append(xp,xp[0]),np.append(yp,yp[0])  #extend arrays
    
    I = 0
    for i in range(N):
        while (I<len(xp)-1):
            if(xp[I]<=x[i]<=xp[I+1] or xp[I+1]<=x[i]<=xp[I]): break
            else: I += 1
        a = (yp[I+1]-yp[I])/(xp[I+1]-xp[I])
        b = yp[I+1]-a*xp[I+1]
        y[i] = a*x[i]+b
    y[N] = y[0]
    
    panel = np.empty(N,dtype=object)
    for i in range(N):
        panel[i] = Panel(x[i],y[i],x[i+1],y[i+1])
        
    return panel

N=20
panel = definePanels(N,xp,yp)
NA = N+1

# plotting the geometry with the panels
valX,valY = 0.2,0.4
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)
plt.plot(np.append([p.xa for p in panel],panel[0].xa),\
        np.append([p.ya for p in panel],panel[0].ya),\
        marker='D',markersize=6,color='r');

A = np.zeros((NA,NA+1),dtype=float)
B = np.zeros_like(A)

CO1 = [p.xc for p in panel]
CO2 = [p.yc for p in panel]

TH = [p.TH for p in panel]
DL = [p.length for p in panel]

HOLDA = 0.0
HOLDB = 0.0

for i in range(N):
    for j in range(N):
        XT = CO1[i] - panel[j].xa
        YT = CO2[i] - panel[j].ya
        X2T = panel[j].xb - panel[j].xa
        Y2T = panel[j].yb - panel[j].ya
                
        X = XT*cos(TH[j]) + YT*sin(TH[j])  # collocation point
        Y = -XT*sin(TH[j]) + YT*cos(TH[j]) # collocation point
        X1 = 0
        Y1 = 0
        X2 = X2T*cos(TH[j]) + Y2T*sin(TH[j])
        Y2 = 0
                
        # Find the length of r1,r2,theta1 and theta2
        R1 = sqrt((X-X1)**2 + (Y-Y1)**2)   # length from panel point 1 to collocation point in panel coords
        R2 = sqrt((X-X2)**2 + (Y-Y2)**2)   # length from panel point 2 to collocation point in panel coords
                
        TH1 = atan2(Y-Y1,X-X1)
        TH2 = atan2(Y-Y2,X-X2)
                
        if i==j:
            Y = 0
            TH1 = 0
        # Compute velocity components as functions of Gamma1 and Gamma2
        # Velocity of panel j due to collocation point i 
        if(i==j):
            U1L = -0.5*(X-X2)/X2*5
            U2L = 0.5*X/X2*5
            W1L = -0.15916*5. 
            W2L = 0.15916 *5.
 #           W1L = -1.
 #           W2L = 1.      
        else:
            U1L = -(Y*log(R2/R1)+X*(TH2-TH1)-X2*(TH2-TH1))/(6.28319*X2)
            U2L = (Y*log(R2/R1) + X*(TH2-TH1))/(6.28319*X2)
            W1L = -((X2-Y*(TH2-TH1)) - X*log(R1/R2) + X2*log(R1/R2))/(6.28319*X2)
            W2L = ((X2 - Y*(TH2-TH1))-X*log(R1/R2))/(6.28319*X2)   
                
        # Transform the local velocities into global velocity function
        U1 = U1L*cos(-TH[j]) + W1L*sin(-TH[j]) 
        U2 = U2L*cos(-TH[j]) + W2L*sin(-TH[j]) 
        W1 = -U1L*sin(-TH[j]) + W1L*cos(-TH[j]) 
        W2 = -U2L*sin(-TH[j]) + W2L*cos(-TH[j]) 
            # Compute the coefficients of gamma in the influence matrix
        if (j==0):
            A[i,0] = -U1*sin(TH[i]) + W1*cos(TH[i])
            HOLDA = -U2*sin(TH[i]) + W2*cos(TH[i])
            B[i,0] = U1*cos(TH[i]) + W1*sin(TH[i])
            HOLDB = U2*cos(TH[i]) + W2*sin(TH[i])
                    
        elif (j==N-1):
            A[i,N-1] = -U1*sin(TH[i]) + W1*cos(TH[i]) + HOLDA
            A[i,N] = -U2*sin(TH[i]) + W2*cos(TH[i])
            B[i,N-1] = U1*cos(TH[i]) + W1*sin(TH[i]) + HOLDB
            B[i,N] = U2*cos(TH[i]) + W2*sin(TH[i])
        else:
            A[i,j] = -U1*sin(TH[i]) + W1*cos(TH[i]) + HOLDA
            HOLDA = -U2*sin(TH[i]) + W2*cos(TH[i])
            B[i,j] = U1*cos(TH[i]) + W1*sin(TH[i]) + HOLDB
            HOLDB = U2*cos(TH[i]) + W2*sin(TH[i])
    
                        
    A[i,NA] = cos(AL) * sin(TH[i]) - sin(AL)* cos(TH[i])
    #A[i,NA] = -freestream.Uinf*cos(freestream.alpha-TH[i])


A[NA-1,0] = 1
A[NA-1,NA-1] = 1

def ToReducedRowEchelonForm(M):
    #if not M : return
    lead = 0
    rowCount = len(M)
    columnCount = len(M[0])
    for r in range(rowCount):
        if lead >= columnCount:
            return
        i=r
        while M[i][lead] == 0:
            i+= 1
            if i==rowCount:
                i = r
                lead += 1
                if columnCount == lead:
                    return
        M[i],M[r] = M[r],M[i]
        lv = M[r][lead]
        M[r] = [mrx/lv for mrx in M[r]]
        for i in range(rowCount):
            if i != r:
                lv = M[i][lead]
                M[i] = [iv - lv*rv for rv,iv in zip(M[r],M[i])]
        lead += 1
   
    return M

R = ToReducedRowEchelonForm(A)
G = R[:,NA]
CL = 0
CP = np.zeros(N,dtype=float)
# calculate variables of interest
for i in range(N):
    VEL = 0
    for j in range(NA):
        VEL = VEL + B[i,j]*G[j]
    
    V = VEL + (cos(AL)*cos(TH[i]) + sin(AL) * sin(TH[i]))
    #V = VEL  -freestream.Uinf*cos(freestream.alpha-TH[i])
    CP[i] = 1 - V**2
    CL = CL -1.0*CP[i]*(cos(AL)*cos(TH[i]) + sin(AL)* sin(TH[i]))*DL[i]




for i in range(N):
    panel[i].Cp = CP[i]
    
print CL
valX,valY = 0.2,0.4
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
Cpmin,Cpmax = min([p.Cp for p in panel]),max([p.Cp for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = Cpmin-valY*(Cpmax-Cpmin),Cpmax+valY*(Cpmax-Cpmin)      
plt.figure(figsize=(10,6))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('$C_p$',fontsize=16)
plt.plot([p.xc for p in panel if p.loc=='extrados'],\
		[p.Cp for p in panel if p.loc=='extrados'],\
		'ro-',linewidth=2)
plt.plot([p.xc for p in panel if p.loc=='intrados'],\
		[p.Cp for p in panel if p.loc=='intrados'],\
		'bo-',linewidth=1)
plt.legend(['extrados','intrados'],'best',prop={'size':14})
#plt.plot([p.xc for p in panel],[p.Cp for p in panel],'ro',linewidth=2)
#plt.plot([p.xc for p in panel],[p.Cp for p in panel],'o-',linewidth=1)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.gca().invert_yaxis();

plt.show();