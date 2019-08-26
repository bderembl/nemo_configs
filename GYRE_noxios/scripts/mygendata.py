#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf

dir0 = "../EXP00/"

file_s = "surface_var.nc"

######## GRID PARAMETERS

si_y = 102
si_x = 102

si_x1 = si_x + 1
si_y1 = si_y + 1

# in m
Lx = 5000.0e3
Ly = 5000.0e3

dx = Lx/si_x;
dy = Ly/si_y;

xx = Lx*(np.arange(0,si_x) + 0.5)/(1.0*si_x)
yy = Ly*(np.arange(0,si_y) + 0.5)/(1.0*si_y)

xx1 = Lx*(np.arange(0,si_x+1) )/(1.0*si_x)
yy1 = Ly*(np.arange(0,si_y+1) )/(1.0*si_y)


xg,yg = np.meshgrid(xx,yy) 
xu,yu = np.meshgrid(xx1[:-1],yy) 
xv,yv = np.meshgrid(xx,yy1[:-1]) 
xc,yc = np.meshgrid(xx1,yy1) 

dx1 = dx*np.ones((si_x))
dy1 = dy*np.ones((si_y))

dz1 = np.array([10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,500,500,500])
np.savetxt(dir0 + 'dz.dat',dz1)

####### forcing and initial conditions

dt = 20
t0 = 2

sst = dt*(1-yg/Ly) + t0


########### create netcdf file

fileout = dir0 + file_s
f = netcdf.netcdf_file(fileout,'w')

f.createDimension('t',None)
f.createDimension('y',si_y)
f.createDimension('x',si_x)

tpo = f.createVariable('t', 'f', ('t',))
ypo = f.createVariable('y', 'f', ('y',))
xpo = f.createVariable('x', 'f', ('x',))

#ssto  = f.createVariable('sst' , 'f', ('t','y','x',))
ssto  = f.createVariable('sst' , 'f', ('y','x',))

ypo[:] = np.arange(si_y)
xpo[:] = np.arange(si_x)

tpo[0] = 0
ssto [:,:] = sst[:,:]

f.close()

