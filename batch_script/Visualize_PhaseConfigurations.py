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
import ctypes
from ctypes import *
import h5py






beta_low=0.3
beta_high=0.35
rho=1
alpha=1
eta1=0
eta2=10
nMAX=30
nbeta=64
e=0
h=1
nu=0
Hinit=1
L=20
MaxP=129 #number of discretized values of the phase

folder_out=("/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_%s/nu_%s/eta2_%s/h_%s" %(e, nu, eta2, h))
beta=np.zeros((nbeta))

num_selected_beta=8
selected_beta=np.arange(0, nbeta, int(nbeta/num_selected_beta))

class Psi_struct(Structure):
	_fields_=[('t1', c_int),
			  ('t2', c_int)]

BASEDIR=("%s/L%s_rho%s_alpha%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s"  %(folder_out, L, rho, alpha, eta1, eta2, e, h, nu, beta_low, beta_high, nMAX, Hinit))

fig1, ax1 = plt.subplots(nrows=8 , ncols=3, sharex=True, sharey=True, figsize=(6, 24) )
#fig1.suptitle(r"$\beta_i=%s$; Initialization:%s" %(b, Hinit))

for b in range(len(selected_beta)):
	WORKDIR=("%s/beta_%s" %(BASEDIR, selected_beta[b]))
#	list_files=fmatch.filter(os.listdir('%s' %WORKDIR), 'Psi_n*.bin' )
#	MC_times=[list_files[d].split('_n') for d in range(len(list_files))]
#	print(MC_times)
	filename_init=("%s/Psi_n0.bin" %WORKDIR)	
	filename_inter=("%s/Psi_n200000.bin" %WORKDIR)
	filename_fin=("%s/Psi_final.bin" %WORKDIR)

	with open('%s' %filename_fin, 'rb') as file:
		result=[]
		data= Psi_struct()
		while file.readinto(data) == sizeof(data):
			result.append((data.t1, data.t2))

	result=np.asarray(result)
	result= result*2*math.pi/MaxP
	phase_field1=result[:,0]
	phase_field2=result[:,1]
	phase_field12=phase_field1 -phase_field2
	phase_field1=phase_field1.reshape(L,L,L)
	phase_field2=phase_field2.reshape(L,L,L)
	phase_field12=phase_field12.reshape(L,L,L)

	x=np.arange(0, L, 1)
	y=np.arange(0, L, 1)
	z=0


	ax1[b,0].set_xticks(np.arange(0, L, step=int(L/4)))
	ax1[b,0].set_yticks(np.arange(0, L, step=int(L/4)))
	if b==0:
		ax1[b,0].set_title('$\phi_1$')
		ax1[b,1].set_title('$\phi_2$')
		ax1[b,2].set_title('$\phi_{12}$')

	im0=ax1[b,0].pcolormesh(x,y, phase_field1[:,:,z], vmin=0, vmax=2*math.pi, cmap='hsv')
	im1=ax1[b,1].pcolormesh(x,y, phase_field2[:,:,z], vmin=0, vmax=2*math.pi, cmap='hsv')
	im2=ax1[b,2].pcolormesh(x,y, phase_field12[:,:,z], vmin=0, vmax=2*math.pi, cmap='hsv')
	cbar=fig1.colorbar(im2, ax=ax1[b,2], ticks=[0, math.pi/2, math.pi, 3*math.pi/2, 2*math.pi] )
	cbar.ax.set_yticklabels([0, '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$' ])
	cbar=fig1.colorbar(im1, ax=ax1[b,1], ticks=[0, math.pi/2, math.pi, 3*math.pi/2, 2*math.pi] )
	cbar.ax.set_yticklabels([0, '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$' ])
	cbar=fig1.colorbar(im0, ax=ax1[b,0], ticks=[0, math.pi/2, math.pi, 3*math.pi/2, 2*math.pi] )
	cbar.ax.set_yticklabels([0, '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$' ])
fig1.tight_layout()
fig1.subplots_adjust(top=0.95, bottom=0.05, hspace=0.2)
plt.show()