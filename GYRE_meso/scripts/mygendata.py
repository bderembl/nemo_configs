#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.netcdf as netcdf
import glob
import spoisson 
from scipy import interpolate

dir0 = "../EXP00/"

file_s = "surface_var.nc"
file_r = "resto.nc"
file_nu = "eddy_diffusivity_2D.nc"

# GRID PARAMETERS
#----------------

si_y = 125 + 2
si_x = 125 + 2
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

# old grid

xxo = Lx*(np.arange(0,N) + 0.5)/(1.0*N)
yyo = Ly*(np.arange(0,N) + 0.5)/(1.0*N)

xx1o = Lx*(np.arange(0,N+1) )/(1.0*N)
yy1o = Ly*(np.arange(0,N+1) )/(1.0*N)

xg,yg = np.meshgrid(xxo,yyo) 
xco,yco = np.meshgrid(xx1o,yy1o) 

b = np.fromfile(allfilesb[-1],'f4').reshape(nl2,N1,N1).transpose(0,2,1)
uv = np.fromfile(allfilesu[-1],'f4').reshape(2*nl2,N1,N1).transpose(0,2,1)


# # interpolate
# # old grid
# Delta = 1/N
# x_old  = np.linspace(0.5*Delta, 1-0.5*Delta,N)
# x2_old = np.linspace(-0.5*Delta, 1+0.5*Delta,N2)

# def bc(psi1,psi2):
  
#   ## to be finished
#   po2 = np.zeros((nl,N2,N2))


#   fr2[:,1:-1,1:-1] = fr[:,1:,1:]
  
#   # BC
#   fr2[:,0,:]  = fr2[:,1,:]
#   fr2[:,-1,:] = fr2[:,-2,:]
#   fr2[:,:,0]  = fr2[:,:,1]
#   fr2[:,:,-1] = fr2[:,:,-2]
  
#   # corners
#   fr2[:,0,0]   = fr2[:,1,1]
#   fr2[:,-1,0]  = fr2[:,-2,1]
#   fr2[:,0,-1]  = fr2[:,1,-2]
#   fr2[:,-1,-1] = fr2[:,-2,-2]


b = Thetas*(b[1:-1,1:,1:] - b.min()) + 2.0
u = Us*uv[2:-2:2,1:,1:]
v = Us*uv[3:-2:2,1:,1:]

ff = fnot + beta*yco[:-1,:-1] # should be at u and v points

# compute pressure for SSH
dudy = np.diff(ff[np.newaxis,:,:]*u,1,1)/dy
dvdx = np.diff(ff[np.newaxis,:,:]*v,1,2)/dx

vort = dvdx[0,:-1,:] - dudy[0,:,:-1]

psi = spoisson.sol(vort[:]) 
psi = psi.reshape((N-1,N-1))
psi = psi*dx*dx/gg # assume dx = dy

psi2 = np.zeros((N, N))
psi2[1:,1:] = psi # not quite right
eta = np.zeros((N+2,N+2))
eta[1:-1,1:-1] = psi2

# interpolate
uvel_n  = np.zeros((si_z,si_y,si_x))
vvel_n  = np.zeros((si_z,si_y,si_x))
theta_n = np.zeros((si_z,si_y,si_x))

eta_n = np.zeros((si_y,si_x))

for nz in range(0,si_z):
  fint = interpolate.interp2d(xxo[:], yyo[:],u[nz,:,:], kind='cubic')
  uvel_n[nz,:,:] = fint(xx,yy)
  
  fint = interpolate.interp2d(xxo, yyo,v[nz,:,:], kind='cubic')
  vvel_n[nz,:,:] = fint(xx,yy)

  fint = interpolate.interp2d(xxo, yyo,b[nz,:,:], kind='cubic')
  theta_n[nz,:,:] = fint(xx,yy)

fint = interpolate.interp2d(xxo, yyo,psi2, kind='cubic')
eta_n = fint(xx,yy)


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
to[:,:,:] = theta_n

uo[:,:,:] = 0.
vo[:,:,:] = 0.
uo[:,:,:] = uvel_n
vo[:,:,:] = vvel_n

eo[:,:] = eta_n

f.close()


# create netcdf file (restoring)
#-------------------------------

tau_dmp = 3600  # restoring time scale

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

ro[:,:,:] = 1/tau_dmp

f.close()


# create netcdf file (viscosity)
#-------------------------------

fileout = dir0 + file_nu
f = netcdf.netcdf_file(fileout,'w')


f.createDimension('y',si_y)
f.createDimension('x',si_x)

ypo = f.createVariable('y', 'd', ('y',))
xpo = f.createVariable('x', 'd', ('x',))

ahmt  = f.createVariable('ahmt_2d' , 'd', ('y','x',))
ahmf  = f.createVariable('ahmf_2d' , 'd', ('y','x',))

ypo[:] = np.arange(si_y)
xpo[:] = np.arange(si_x)

ahmt[:,:] = 100
ahmf[:,:] = 100

f.close()
