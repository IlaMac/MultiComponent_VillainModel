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

beta_low=0.42
beta_high=0.44
rho=1
alpha=1
eta1=0
eta2=10
nMAX=30
nbeta=64
e=0.
h=1
nu=0.4
Hinit=1
L=[32]

if( (nu).is_integer()): nu=int(nu)
if( (e).is_integer()): e=int(e)

folder_out=("/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_%s/nu_%s/eta2_%s/h_%s" %(e, nu, eta2, h))
beta=np.zeros((nbeta))


for l in range(len(L)):
	N=L[l]**3

	BASEDIR=("%s/L%s_rho%s_alpha%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s"  %(folder_out, L[l], rho, alpha, eta1, eta2, e, h, nu, beta_low, beta_high, nMAX, Hinit))
	folder_ranks=("%s/Ranks" %BASEDIR)
	try:
	    os.makedirs(folder_ranks)    
	    print("Directory " , folder_ranks ,  " Created ")
	except FileExistsError:
	    print("Directory " , folder_ranks ,  " already exists") 
	temp_length=[]
	for b in range(nbeta):
		file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
		b_rank=np.asarray(file['Measurements']['rank'])
		temp_length.append(len(b_rank))
	print(temp_length, min(temp_length))
	min_length= min(temp_length)
	ranks=np.zeros((nbeta, min_length))

	beta_rank=[]
	for b in range(nbeta):
		beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)

		file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
		b_rank=np.asarray(file['Measurements']['rank'])

		for r in range(min_length):
			#print(b, r, np.shape(ranks), ranks)
			ranks[b_rank[r], r]=beta[b]
#

	for b in range(nbeta):
		print(b, ranks[b, :])
		f, (ax1) = plt.subplots(ncols=1, figsize=(6,6) )
		f.suptitle(r"$rank=%s$" %(b))
		ax1.set_xlabel(r"$t_{MC}$")
		ax1.set_ylim((beta[0], beta[nbeta-1]))
		ax1.plot(ranks[b, :])
		f.savefig("%s/rank_%s.png" %(folder_ranks, b))
		plt.show()
#




