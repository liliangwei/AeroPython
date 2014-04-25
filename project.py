# Project 
# Starson
# Two element panel method

import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt

alpha = 5.       # the angle of attack from main elemnt chord line
Uinf = 20.

AL = alpha / 57.2958     # get into radians


# read of the geometry
#coords = np.loadtxt(fname='C:/Users/llwei89/Documents/Github/AeroPython/resources/n0012.dat')
coords = np.loadtxt(fname='/home/starson/AeroPython/resources/s1223.dat')
xpM,ypM = coords[:,0],coords[:,1]            # read in the original airfoil

# Scale the second flap
xpF = 0.5*np.copy(xpM)                       # scaled x coordinate of flap
ypF = 0.5*np.copy(ypM)                       # scaled y coordinate of flap

# rotate the flap angle delta 35 degree down
delta = 35

x1 = np.zeros_like(xpF)      
y1 = np.zeros_like(ypF)

L = len(xpM)

for i in range(L):
    r  = sqrt(xpF[i]**2+ypF[i]**2)
    thetatr = np.arctan2(ypF[i],xpF[i])
    theta = thetatr - delta*pi/180.0
    x1[i] = r*cos(theta)
    y1[i] = r*sin(theta)

xpF = x1
ypF = y1

xpF = xpF + 0.92                            # translated x coordinate
ypF = ypF -0.06                             # translated y coordinate

# plotting the geometry
valX,valY = 0.2,0.4                          #value for plot margin    
xmin,xmax = min(xpM),max(xpF)
ymin,ymax = min(ypF),max(ypM)
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)

size=10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xpM,ypM,xpF,ypF,'k-',linewidth=2);

# Create panels for both main and flap airfoils.

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
    

NM,NF = 100,50                                     # Number of panels for each


Ntotal = NM + NF                                  # Total number of panels
NA = Ntotal + 2                                   # Size of A matrix

panelM = definePanels(NM,xpM,ypM)
panelF = definePanels(NF,xpF,ypF)

# plotting the geometry with the panels
valX,valY = 0.2,0.4
xmin,xmax = min([p.xa for p in panelM]),max([p.xa for p in panelF])
ymin,ymax = min([p.ya for p in panelF]),max([p.ya for p in panelM])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xpM,ypM,xpF,ypF,'k-',linewidth=2)
plt.plot(np.append([p.xa for p in panelM],panelM[0].xa),\
        np.append([p.ya for p in panelM],panelM[0].ya),\
        np.append([p.xa for p in panelF],panelF[0].xa),\
        np.append([p.ya for p in panelF],panelF[0].ya),\
        marker='D',markersize=6,color='r');
        
# Create A matrix
A = np.zeros((NA,NA+1),dtype=float)
B = np.zeros_like(A)
# Merge all coefficients into one set
# Do flap then main

CO1 = np.append([p.xc for p in panelF],[p.xc for p in panelM])
CO2 = np.append([p.yc for p in panelF],[p.yc for p in panelM])

TH = np.append([p.TH for p in panelF],[p.TH for p in panelM])
DL = np.append([p.length for p in panelF],[p.length for p in panelM])


for i in range(Ntotal):

#determine if we are on the flap or main
    if(i<NF):
    # Dealing with the both collocation points on flap

        for j in range(Ntotal+1):
            if(j<NF):
        # We are dealing with both collocation points on flap
        # We find the influence coefficient for a specific collocation 
        # point, i, on each of the panels, j. Then we move to the next
        # collocation point

        # Convert collocation point to local panel coordinates   
                XT = CO1[i] - panelF[j].xa
                YT = CO2[i] - panelF[j].ya
                X2T = panelF[j].xb - panelF[j].xa
                Y2T = panelF[j].yb - panelF[j].ya
                
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
                    U1L = -0.5*(X-X2)/X2
                    U2L = 0.5*X/X2
                    W1L = -Uinf
                    W2L = Uinf
                    
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
                if ((j==0) or (j==NF)):
                    A[i,0] = -U1*sin(TH[i]) + W1*cos(TH[i])
                    HOLDA = -U2*sin(TH[i]) + W2*cos(TH[i])
                    B[i,0] = U1*cos(TH[i]) + W1*sin(TH[i])
                    HOLDB = U2*cos(TH[i]) + W2*sin(TH[i])
                    
                elif (j==NF-1):
                    A[i,NF-1] = -U1*sin(TH[i]) + W1*cos(TH[i]) + HOLDA
                    A[i,NF] = -U2*sin(TH[i]) + W2*cos(TH[i])
                    B[i,NF-1] = U1*cos(TH[i]) + W1*sin(TH[i]) + HOLDB
                    B[i,NF] = U2*cos(TH[i]) + W2*sin(TH[i])
                else:
                    A[i,j] = -U1*sin(TH[i]) + W1*cos(TH[i]) + HOLDA
                    HOLDA = -U2*sin(TH[i]) + W2*cos(TH[i])
                    B[i,j] = U1*cos(TH[i]) + W1*sin(TH[i]) + HOLDB
                    HOLDB = U2*cos(TH[i]) + W2*sin(TH[i])
            elif (j>NF):
        # collocation point i is on the flap, collocation point j is on main
        
        # we find the influence coefficient for a specific collocation
        # point, i, on each of panels, j. Then we move on to the next collocation point
        
        # since there are two separate airfoils, point j=NF+1
        # should not be solved for (this is the 2nd edge of the last panel)
        # During the for loop, any j after this should be j-1, 
        # the reference the correct panel, except when referring to matrix A
        
        # Convert collocation point to local coords
                XT = CO1[i] - panelM[j-NF-1].xa
                YT = CO2[i] - panelM[j-NF-1].ya
                X2T = panelM[j-NF-1].xb-panelM[j-NF-1].xa
                Y2T = panelM[j-NF-1].yb-panelM[j-NF-1].ya
               
                X = XT*cos(TH[j-1]) + YT*sin(TH[j-1])   # collocation point
                Y = -XT*sin(TH[j-1]) + YT*cos(TH[j-1])  # collocation point
                X1 = 0
                Y1 = 0
                X2 = X2T*cos(TH[j-1]) + Y2T*sin(TH[j-1])
                Y2 = 0
        # Find the length of r1,r2, theta1 and theta2
                R1 = sqrt((X-X1)**2 + (Y-Y1)**2)
                R2 = sqrt((X-X2)**2 + (Y-Y2)**2)
                
                TH1 = atan2(Y-Y1,X-X1)
                TH2 = atan2(Y-Y2,X-X2)
                
                if (i==j-1):
                    Y = 0
                    TH1 = 0
                    
                
        # Compute velocity components as functions of Gamma1 and Gamm2
        # velocity of panel j due to collocation point i
                if (i==j-1):
                    U1L = -0.5*(X-X2)/X2
                    U2L = 0.5*X/X2
                    W1L = -Uinf
                    W2L = Uinf
                else:
                    U1L = -(Y*log(R2/R1)+X*(TH2-TH1)-X2*(TH2-TH1))/(6.28319*X2)
                    U2L = (Y*log(R2/R1) + X*(TH2-TH1))/(6.28319*X2)
                    W1L = -((X2-Y*(TH2-TH1))-X*log(R1/R2) + X2*log(R1/R2))/(6.28310*X2)
                    W2L = ((X2-Y*(TH2-TH1))-X*log(R1/R2))/(6.28319*X2)
                    
        # Transferom the local velocities into global velocity function
                U1 = U1L*cos(-TH[j-1]) + W1L*sin(-TH[j-1])
                U2 = U2L*cos(-TH[j-1]) + W2L*sin(-TH[j-1])
                W1 = -U1L*sin(-TH[j-1]) + W1L*cos(-TH[j-1])
                W2 = -U2L*sin(-TH[j-1]) + W2L*cos(-TH[j-1])
        
        # Compute the coefficients of gamma in the influence matrix
                if (j==NF+1):
                    A[i,NF+1] = -U1*sin(TH[i]) + W1*cos(TH[i])
                    HOLDA = -U2*sin(TH[i]) + W2*cos(TH[i])
                    B[i,NF+1] = U1*cos(TH[i]) + W1*sin(TH[i])
                    HOLDB = U2*cos(TH[i]) + W2*sin(TH[i])
                elif (j==Ntotal):
                    A[i,Ntotal] = -U1*sin(TH[i]) + W1*cos(TH[i]) + HOLDA
                    A[i,NA-1] = -U2*sin(TH[i]) + W2*cos(TH[i])
                    B[i,Ntotal] = U1*cos(TH[i]) + W1*sin(TH[i]) + HOLDB
                    B[i,NA-1] = U2*cos(TH[i]) + W2*sin(TH[i])
                else:
                    A[i,j] = -U1*sin(TH[i]) + W1*cos(TH[i]) + HOLDA
                    HOLDA = -U2*sin(TH[i]) + W2*cos(TH[i])
                    B[i,j] = U1*cos(TH[i]) + W1*sin(TH[i]) + HOLDB
                    HOLDB = U2*cos(TH[i]) + W2*sin(TH[i])
         
        A[i,NA] = cos(AL) * sin(TH[i]) -sin(AL)*cos(TH[i])

    else:
    # We are dealing with the collocation point on the main
        for j in range(Ntotal+1):
            if (j<NF):
        # CO point i is on main, CO point j is on flap
        
        # We find the influence coefficient for a specific collocation 
        # point, i, on each of panels, j. Then we move to the n collocation point.
        
        # Convert collocation point to local panel coordinates
                XT = CO1[i] - panelF[j].xa
                YT = CO2[i] - panelF[j].ya
                X2T = panelF[j].xb - panelF[j].xa
                Y2T = panelF[j].yb - panelF[j].ya
                
                X = XT*cos(TH[j]) + YT*sin(TH[j])     # collocation point
                Y = -XT*sin(TH[j]) + YT*cos(TH[j])    # collocation point
                X1 = 0
                Y1 = 0
                X2 = X2T*cos(TH[j]) + Y2T*sin(TH[j])
                Y2 = 0
                
        # find the length of r1, r2, theta1 and theta2
                R1 = sqrt((X-X1)**2 + (Y-Y1)**2)      # length from panel point 1 collocation point , panel coords
                R2 = sqrt((X-X2)**2 + (Y-Y2)**2)      # length from panel point 2 collocation point, in panel coords
                
                TH1 = atan2(Y-Y1,X-X1)
                TH2 = atan2(Y-Y2,X-X2)
                
                if (i==j):
                    Y = 0
                    TH1 = 0
                
        # Compute velocity componets as function of Gamma1 and Gamma2
        # Velocity of panel j due to collocation point 1
                if(i==j):
                    U1L = -0.5*(X-X2)/X2
                    U2L = 0.5*X/X2
                    W1L = -Uinf
                    W2L = Uinf
                else:
                    U1L = -(Y*log(R2/R1)+X*(TH2-TH1)-X2*(TH2-TH1))/(6.28319*X2)
                    U2L = (Y*log(R2/R1) + X*(TH2-TH1))/(6.28319*X2)
                    W1L = -((X2-Y*(TH2-TH1)) - X*log(R1/R2) + X2*log(R1/R2))/(6.28319*X2)
                    W2L = ((X2 - Y*(TH2-TH1))-X*log(R1/R2))/(6.28319*X2)
        # Transform the local velocities into global velocity functions
                
                U1 = U1L*cos(-TH[j]) + W1L*sin(-TH[j])
                U2 = U2L*cos(-TH[j]) + W2L*sin(-TH[j])
                W1 = -U1L*sin(-TH[j]) + W1L*cos(-TH[j])
                W2 = -U2L*sin(-TH[j]) + W2L*cos(-TH[j])
                
                if (j==0):
                    A[i,0] = -U1*sin(TH[i]) + W1*cos(TH[i])
                    HOLDA = -U2*sin(TH[i]) + W2*cos(TH[i])
                    B[i,0] = U1*cos(TH[i]) + W1*sin(TH[i])
                    HOLDB = U2*cos(TH[i]) + W2*sin(TH[i])
                elif(j==NF-1):
                    A[i,NF-1] = -U1*sin(TH[i]) + W1*cos(TH[i]) + HOLDA
                    A[i,NF] = -U2*sin(TH[i]) + W2*cos(TH[i])
                    B[i,NF-1] = U1*cos(TH[i]) + W1*sin(TH[i]) + HOLDB
                    B[i,NF] = U2*cos(TH[i]) * W2*sin(TH[i])
                else:
                    A[i,j] = -U1*sin(TH[i]) + W1*cos(TH[i]) + HOLDA
                    HOLDA = -U2*sin(TH[i]) + W2*cos(TH[i])
                    B[i,j] = U1*cos(TH[i]) + W1*sin(TH[i]) + HOLDB
                    HOLDB = U2*cos(TH[i]) + W2*sin(TH[i])
            
            elif (j>NF):
            
            # collocation point i is on main, collocation point j is on main
            
            # we find the influence coefficient for a specific collocation
            # point, i, on each of panel, j. Then we move on to the next collocation point
            
            # since there are two separate airfolds, point j = NF
            # should not be solved for (this is the 2nd edge of the last panel)
            # During the for loop, any j after
            # this should be j-1, to reference the correct panel, except when referring to matrix A
            
            # Convert collocation point To local panel coordinates
                XT = CO1[i] - panelM[j-NF-1].xa
                YT = CO2[i] - panelM[j-NF-1].ya
                X2T = panelM[j-NF-1].xb - panelM[j-NF-1].xa
                Y2T = panelM[j-NF-1].yb - panelM[j-NF-1].ya
                
                X = XT*cos(TH[j-1]) + YT*sin(TH[j-1])
                Y = -XT*sin(TH[j-1]) + YT*cos(TH[j-1])
                X1 = 0
                Y1 = 0
                X2 = X2T*cos(TH[j-1]) + Y2T*sin(TH[j-1])
                Y2 = 0
                
            # Find the length of r1, r2, theta1, and theta2
                R1 = sqrt((X-X1)**2 + (Y-Y1)**2)  # length from panel point 1 to collocation point, in panel coords
                R2 = sqrt((X-X2)**2 + (Y-Y2)**2)  # length from panel point 2 to collocation point, in panel coords
                
                TH1 = atan2(Y-Y1,X-X1)
                TH2 = atan2(Y-Y2,X-X2)
                
                if i==j-1:
                    Y = 0
                    TH1 = 0
                
            # Compute velocity components as functions of Gamma1 and Gamma2
            # Velocity of panel j due to collocation point i 
                if i==j-1:
                    U1L = -0.5*(X-X2)/X2
                    U2L = 0.5*X/X2
                    W1L = -Uinf
                    W2L = Uinf
                else:
                    U1L = -(Y*log(R2/R1) + X*(TH2-TH1) - X2*(TH2-TH1))/(6.28319*X2)
                    U2L = (Y*log(R2/R1) + X*(TH2-TH1))/(6.28319*X2)
                    W1L = -((X2-Y*(TH2-TH1)) - X*log(R1/R2) + X2*log(R1/R2))/(6.28319*X2)
                    W2L = ((X2-Y*(TH2-TH1)) - X*log(R1/R2))/(6.28319*X2)
                
            # Transform the local velocities into global velocity functions
                U1 = U1L*cos(-TH[j-1]) + W1L*sin(-TH[j-1])
                U2 = U2L*cos(-TH[j-1]) + W2L*sin(-TH[j-1])
                W1 = -U1L*sin(-TH[j-1]) + W1L*cos(-TH[j-1])
                W2 = -U2L*sin(-TH[j-1]) + W2L*cos(-TH[j-1])
                
            # Compute the coefficients of gamma in the influence matrix
                if (j == NF+1):
                    A[i,NF+1] = -U1*sin(TH[i]) + W1*cos(TH[i])
                    HOLDA = -U2*sin(TH[i]) + W2*cos(TH[i])
                    B[i,NF+1] = U1*cos(TH[i]) + W1*sin(TH[i])
                    HOLDB = U2*cos(TH[i]) + W2*sin(TH[i])
                elif (j==Ntotal):
                    A[i,Ntotal] = -U1*sin(TH[i]) + W1*cos(TH[i]) + HOLDA
                    A[i,NA-1] = -U2*sin(TH[i]) + W2*cos(TH[i])
                    B[i,Ntotal] = U1*cos(TH[i]) + W1*sin(TH[i]) + HOLDB
                    B[i,NA-1] = U2*cos(TH[i]) + W2*sin(TH[i])
                else:
                    A[i,j] = -U1*sin(TH[i]) + W1*cos(TH[i]) + HOLDA
                    HOLDA = -U2*sin(TH[i]) + W2*cos(TH[i])
                    B[i,j] = U1*cos(TH[i]) + W1*sin(TH[i]) + HOLDB
                    HOLDB = U2*cos(TH[i]) + W2*sin(TH[i])
                
        A[i,NA] = cos(AL) * sin(TH[i]) - sin(AL)* cos(TH[i])


        # Add both kutta conditions. Be careful of where the ones are
        # matrix columns NF and Ntotal+1 are the last edges of the airfold

#Flap
A[NA-2,0] = 1
A[NA-2,NF] = 1
#Main
A[NA-1,NF+1] = 1
A[NA-1,NA-1] = 1

# Linear solver
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
CP = np.zeros(Ntotal,dtype=float)
# calculate variables of interest
for i in range(Ntotal):
    VEL = 0
    for j in range(NA):
        VEL = VEL + B[i,j]*G[j]
    
    V = VEL + cos(AL)*cos(TH[i]) + sin(AL) * sin(TH[i])
    CP[i] = 1 - V**2
    CL = CL -1.0*CP[i]*(cos(AL)*cos(TH[i]) + sin(AL)* sin(TH[i]))*DL[i]

CP = np.transpose(CP)

for i in range(Ntotal):
    if (i<NF):
        panelF[i].Cp = CP[i]
    else:
        panelM[i-NF].Cp = CP[i]

    CL = CL + -1.0*CP[i]*(cos(AL)*cos(TH[i]) + sin(AL)* sin(TH[i]))*DL[i]


# Plot pressure coefficients
valX,valY = 0.2,0.4
xmin,xmax = min([p.xa for p in panelM]),max([p.xa for p in panelF])
Cpmin,Cpmax = min([p.Cp for p in panelM]),max([p.Cp for p in panelM])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = Cpmin-valY*(Cpmax-Cpmin),Cpmax+valY*(Cpmax-Cpmin)      
plt.figure(figsize=(10,6))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('$C_p$',fontsize=16)
plt.plot([p.xc for p in panelM if p.loc=='extrados'],\
		[p.Cp for p in panelM if p.loc=='extrados'],\
		'ro-',linewidth=2)
plt.plot([p.xc for p in panelM if p.loc=='intrados'],\
		[p.Cp for p in panelM if p.loc=='intrados'],\
		'bo-',linewidth=1)
plt.plot([p.xc for p in panelF if p.loc=='extrados'],\
		[p.Cp for p in panelF if p.loc=='extrados'],\
		'ro-',linewidth=2)
plt.plot([p.xc for p in panelF if p.loc=='intrados'],\
		[p.Cp for p in panelF if p.loc=='intrados'],\
		'bo-',linewidth=1)
plt.legend(['extrados','intrados'],'best',prop={'size':14})
plt.plot([p.xc for p in panelF],[p.Cp for p in panelF],'ro',linewidth=2)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.gca().invert_yaxis();



plt.show()