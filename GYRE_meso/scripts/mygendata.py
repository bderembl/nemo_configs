#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf
import glob
import spoisson 

dir0 = "../EXP00/"

file_s = "surface_var.nc"
file_r = "resto.nc"

# GRID PARAMETERS
#----------------

si_y = 64 + 2
si_x = 64 + 2
si_z = 30

si_x1 = si_x + 1
si_y1 = si_y + 1

# in m
Lx = 5000.0e3
Ly = 5000.0e3
Lz = 5000.

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

dz1 = Lz/si_z*np.ones(si_z)
np.savetxt(dir0 + 'dz.dat',dz1)

# initial conditions
#-------------------

# PG scales
L = 5000e3  # m
H = 5000    # m
beta = 2.0e-11 # 1/m/s
N2 = 1e-6  #  (1/s**2)
Bs = N2*H
Thetas = Bs/10/2e-4 # 1/g alpha
Us = N2*H**2/(beta*L**2)
fnot = 3e-5
gg = 9.80665 # nemo value

ff = fnot + beta*yg # should be at u and v points

dir_data = "data/"

fileb = 'b*'
fileu = 'u*'

allfilesb = sorted(glob.glob(dir_data + fileb));
allfilesu = sorted(glob.glob(dir_data + fileu));
nb_files  = len(allfilesb);

# dimensions
b = np.fromfile(allfilesb[0],'f4')
N = int(b[0])
N1 = N + 1
nl2 = int(len(b)/N1**2)
nl = nl2 - 2

b = np.fromfile(allfilesb[-1],'f4').reshape(nl2,N1,N1).transpose(0,2,1)
uv = np.fromfile(allfilesu[-1],'f4').reshape(2*nl2,N1,N1).transpose(0,2,1)

b = Thetas*(b[1:-1,1:,1:] - b.min()) + 2.0
u = Us*uv[2:-2:2,1:,1:]
v = Us*uv[3:-2:2,1:,1:]


# compute pressure for SSH
dudy = np.diff(ff[np.newaxis,1:-1,1:-1]*u,1,1)/dy
dvdx = np.diff(ff[np.newaxis,1:-1,1:-1]*v,1,2)/dx

vort = dvdx[0,:-1,:] - dudy[0,:,:-1]

psi = spoisson.sol(vort[:]) 
psi = psi.reshape((si_y-3,si_x-3))
psi = psi*dx*dx/gg # assume dx = dy

psi2 = np.zeros((si_y-2, si_x-2))
psi2[1:,1:] = psi # not quite right
eta = np.zeros((si_y,si_x))
eta[1:-1,1:-1] = psi2

# create netcdf file (initial state)
fileout = dir0 + 'istate.nc'
f = netcdf.netcdf_file(fileout,'w')

f.createDimension('z',si_z)
f.createDimension('y',si_y)
f.createDimension('x',si_x)

zpo = f.createVariable('z', 'd', ('z',))
ypo = f.createVariable('y', 'd', ('y',))
xpo = f.createVariable('x', 'd', ('x',))

to  = f.createVariable('votemper' , 'd', ('z','y','x',))
so  = f.createVariable('vosaline' , 'd', ('z','y','x',))
uo  = f.createVariable('vozocrtx' , 'd', ('z','y','x',))
vo  = f.createVariable('vomecrty' , 'd', ('z','y','x',))

eo  = f.createVariable('sossheig' , 'd', ('y','x',))

zpo[:] = np.arange(si_z)
ypo[:] = np.arange(si_y)
xpo[:] = np.arange(si_x)
 
so[:,:,:] = 35.

to[:,:,:] = 0.
to[:,1:-1,1:-1] = b

uo[:,:,:] = 0.
vo[:,:,:] = 0.
uo[:,1:-1,1:-1] = u
vo[:,1:-1,1:-1] = v

eo[:,:] = eta

f.close()


# forcing 
#--------

dt = 20
t0 = 2

sst = dt*(1-yg/Ly) + t0


# create netcdf file (surface forcing)
fileout = dir0 + file_s
f = netcdf.netcdf_file(fileout,'w')

f.createDimension('y',si_y)
f.createDimension('x',si_x)

ypo = f.createVariable('y', 'f', ('y',))
xpo = f.createVariable('x', 'f', ('x',))

#ssto  = f.createVariable('sst' , 'f', ('t','y','x',))
ssto  = f.createVariable('sst' , 'f', ('y','x',))

ypo[:] = np.arange(si_y)
xpo[:] = np.arange(si_x)

ssto [:,:] = sst[:,:]

f.close()


# create netcdf file (restoring)
#-------------------------------

fileout = dir0 + file_r
f = netcdf.netcdf_file(fileout,'w')


f.createDimension('z',si_z)
f.createDimension('y',si_y)
f.createDimension('x',si_x)

zpo = f.createVariable('z', 'd', ('z',))
ypo = f.createVariable('y', 'd', ('y',))
xpo = f.createVariable('x', 'd', ('x',))

ro  = f.createVariable('resto' , 'd', ('z','y','x',))

zpo[:] = np.arange(si_z)
ypo[:] = np.arange(si_y)
xpo[:] = np.arange(si_x)

ro[:,:,:] = 1e-3

f.close()
