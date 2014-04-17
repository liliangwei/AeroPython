# Project 
# Starson
# Two element panel method

import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt

# read of the geometry
#coords = np.loadtxt(fname='C:/Users/llwei89/Documents/Github/AeroPython/resources/n0012.dat')
coords = np.loadtxt(fname='/home/starson/AeroPython/resources/n0012.dat')
xp1,yp1 = coords[:,0],coords[:,1]

# Scale the second flap
xp2 = 0.5*xp1
yp2 = 0.5*yp2

# rotate the flap angle delta 35 degree down
delta = 35

x1,y1 = np.zeros(xp2,yp2)
for i in range(len(xp2)):
    r  = np.sqrt(xp2(i)**2+yp2(i)**2)
    thetatr = np.arctan2(yp2