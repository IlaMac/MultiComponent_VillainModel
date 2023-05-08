import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys
import os
import math
from statsmodels.graphics.tsaplots import plot_acf
import statsmodels.api as sm
from statsmodels.tsa.stattools import acf
import scipy.integrate as integrate
import random
import h5py
from scipy.optimize import fsolve
from scipy.optimize import bisect
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

plt.rc('text',usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{bm}')
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.serif'] = 'Computer Modern'
plt.rcParams['axes.linewidth']  = 3.0
plt.rcParams['axes.labelsize']  = 30
plt.rcParams.update({'font.size': 22})
plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25
plt.rcParams['xtick.major.size'] = 4
plt.rcParams['ytick.major.size'] = 4
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['legend.fontsize']  = 22
plt.rcParams['legend.frameon']  = True
plt.rcParams["legend.fancybox"] = True
plt.rcParams["legend.shadow"] = True
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["axes.facecolor"] = 'white'
plt.rcParams["axes.edgecolor"] = 'black'

plt.rcParams['figure.figsize'] = 8, 8

plt.rcParams['xtick.major.width'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['xtick.minor.width'] = 0
plt.rcParams['ytick.minor.width'] = 0



bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)


import ctypes
from ctypes import *

class YourStruct(Structure):
    _fields_ = [('t', c_int)]
    
    
    
################################################

rho=1
eta1=0
eta2=0
e=0
h=1
nu=-0.75
beta_low=0.1
beta_high=0.5
nMAX=30
H_init=0
MAXP=129
nbeta=64 

dp=2*math.pi/MAXP
###############################################


color=iter(plt.cm.rainbow(np.linspace(0,1,6)))

HOMEDIR=("/Users/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/e_%s/nu_%s/h_%s" %(e, nu, h))
folder_Fig=("/Users/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/e_%s/nu_%s/h_%s" %(e, nu, h))


L=8
beta= np.zeros((nbeta))
for b in range(nbeta):
  beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)

  BASEDIR=("%s/L%d_rho%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s/beta_%s/"  %(HOMEDIR, L, rho, eta1, eta2, e,  h, nu, beta_low, beta_high, nMAX, H_init, b))
  with open('%s/Psi_n100000.bin' %(BASEDIR), 'rb') as file:
      result = []
      data = YourStruct()
      while file.readinto(data) == sizeof(data):
          result.append((data.t))
  
  print(data.t)
  
  result=np.asarray(result)
  print(np.shape(result), len(result))
  theta_field= dp*result[:]
  print(np.shape(theta_field), len(theta_field))
  theta_field=theta_field.reshape(L*L*L,2)
  theta_field_0=theta_field[:,0]
  theta_field_0=theta_field_0.reshape(L,L,L)
  theta_field_1=theta_field[:,1]
  theta_field_1=theta_field_1.reshape(L,L,L)
  print(np.shape(theta_field_0))
  print(theta_field_0[0,0,0])
  

  
  fig = plt.figure(figsize=(15,15))
  ax = fig.gca(projection='3d')
  #ax.set_title(r"Component 1")
  
  x, y, z = np.meshgrid(np.arange(0, L, 1),
                        np.arange(0, L, 1),
                        np.arange(0, L, 1))
  
  #theta_field_0=theta_field_0[:,:,:]*0+0.5*math.pi
  
  u = np.cos(theta_field_0[x,y,z]) 
  v = np.sin(theta_field_0[x,y,z]) 
  w = (theta_field_0[x,y,z]*0)
  
  cmap='hsv'
  # Color by azimuthal angle
  c=0
  c = np.arctan2(v, u)+ math.pi
  # Flatten and normalize
  c=c.ravel()
  print(c.mean())
  
  # Repeat for each body line and two head lines
  c = np.concatenate((c, np.repeat(c, 2)))
  norm = mcolors.Normalize(vmin=0, vmax=2*math.pi,clip=False)
  # Colormap
  c = getattr(plt.cm, cmap)(norm(c))
  
  ax.annotate("Component 1" + "\n" + r"$\beta=%s;$ $e=%s$; $\eta_1=%s$; $\eta_2=%s$; $\nu=%s$" %(beta[b], e, eta1, eta2, nu) , (0.5, 0.97), xycoords='axes fraction',  ha="center", va="center", bbox=bbox_props)
  #ax.annotate("Component 1" + "\n" + r"$\beta=%s$" %(beta_high)+ "\n" + r"$e=%s$" %e + "\n" + r"$\eta=%s$"  %eta  + "\n" + r"$\nu=%s$" %nu  , (0.15, 0.85), xycoords='axes fraction',  ha="center", va="center", bbox=bbox_props)
  
  q = ax.quiver(x, y, z, u, v, w, colors=c, cmap=cmap,  lw=2, length=0.5)
                #length=0.5,  norm=norm, lw=2, colors=c)
  q.set_edgecolor(c)
  q.set_facecolor(c)
  cbar=fig.colorbar(q,fraction=0.036, pad=0.04, ticks=[0, 0.25,0.5, 0.75,1 ])
  cbar.ax.set_yticklabels(['0', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'])  # vertically oriented colorbar
  
  fig.savefig("%s/SpinPhoto_L%s_c1_e%s_eta1%s_eta2%s_nu%s_beta%s.png" %(folder_Fig, L, e, eta1, eta2, nu, b))
  
  fig1 = plt.figure(figsize=(15,15))
  ax1 = fig1.gca(projection='3d')
  #ax.set_title(r"Component 1")
  
  x, y, z = np.meshgrid(np.arange(0, 8, 1),
                        np.arange(0, 8, 1),
                        np.arange(0, 8, 1))
  
  #theta_field_0=theta_field_0[:,:,:]*0+0.5*math.pi
  
  u = np.cos(theta_field_1[x,y,z]) 
  v = np.sin(theta_field_1[x,y,z]) 
  w = (theta_field_1[x,y,z]*0)
  
  cmap='hsv'
  # Color by azimuthal angle
  c=0
  c = np.arctan2(v, u)+ math.pi
  # Flatten and normalize
  c=c.ravel()
  print(c.mean())
  
  # Repeat for each body line and two head lines
  c = np.concatenate((c, np.repeat(c, 2)))
  norm = mcolors.Normalize(vmin=0, vmax=2*math.pi,clip=False)
  # Colormap
  c = getattr(plt.cm, cmap)(norm(c))
  
  ax1.annotate("Component 2" + "\n" + r"$\beta=%s;$ $e=%s$; $\eta_1=%s$; $\eta_2=%s$; $\nu=%s$" %(beta[b], e, eta1, eta2,  nu) , (0.5, 0.97), xycoords='axes fraction',  ha="center", va="center", bbox=bbox_props)
  #ax.annotate("Component 1" + "\n" + r"$\beta=%s$" %(beta_high)+ "\n" + r"$e=%s$" %e + "\n" + r"$\eta=%s$"  %eta  + "\n" + r"$\nu=%s$" %nu  , (0.15, 0.85), xycoords='axes fraction',  ha="center", va="center", bbox=bbox_props)
  
  q = ax1.quiver(x, y, z, u, v, w, colors=c, cmap=cmap,  lw=2, length=0.5)
                #length=0.5,  norm=norm, lw=2, colors=c)
  q.set_edgecolor(c)
  q.set_facecolor(c)
  cbar1=fig1.colorbar(q,fraction=0.036, pad=0.04, ticks=[0, 0.25,0.5, 0.75,1 ])
  cbar1.ax.set_yticklabels(['0', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'])  # vertically oriented colorbar
  
  fig1.savefig("%s/SpinPhoto_L%s_c2_e%s_eta1%s_eta2%s_nu%s_beta%s.png" %(folder_Fig, L, e, eta1, eta2, nu, b))
  
