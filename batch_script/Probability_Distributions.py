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

L=8
rho=1

e=0
nu=0
eta=0
a=h=1
H_init=0
beta_low=0.1
beta_high=0.5
beta_tgt=0.45
nrw=0
nbl=20

nbetas=64
beta=np.zeros((nbetas))

Nbins=50
transient=0.2

folder_out=("/Users/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/e_%s/nu_%s/h_%s" %(e, nu, a))
BASEDIR=("%s/L%d_rho%s_eta%s_e%s_h%s_nu%s_bmin%s_bmax%s_init%s"  %(folder_out, L, rho, eta, e,  a, nu, beta_low, beta_high, H_init))


U_mean=np.zeros(nbetas)
U_var=np.zeros(nbetas)

M_mean=np.zeros(nbetas)
M_var=np.zeros(nbetas)

E_mean=np.zeros(nbetas)
E_var=np.zeros(nbetas)



for b in range(nbetas):
	beta[b]=beta_low +b*(beta_high -beta_low)/(nbetas-1)

	file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
	M=np.asarray(file['Measurements']['m_phase'])
	print(len(M[:,0]))
	M=M[:,0]
	M=M[int(transient*len(M)):]
	M_mean[b]=np.mean(M)
	M_var[b]=np.var(M)

	U=np.asarray(file['Measurements']['U'])
	U=U[int(transient*len(U)):]

	E=np.asarray(file['Measurements']['E'])
	E=beta[b]*E/((L*a)**3)
	E=E[int(transient*len(E)):]

	histM, bin_edgesM = np.histogram(M)
	PM = sm.nonparametric.KDEUnivariate(M)
	PM.fit()

	histU, bin_edgesU = np.histogram(U)
	PU = sm.nonparametric.KDEUnivariate(U)
	PU.fit()

	histE, bin_edgesE = np.histogram(E)
	PE = sm.nonparametric.KDEUnivariate(E)
	PE.fit()

	#plt.rc('text', usetex=True)
	#plt.rc('font', family='serif')




	#plt.rc('text.latex', preamble=r'\usepackage{bm}')
	#fig, (axs) = plt.subplots(nrows=3, ncols=3, figsize=(8,8))
	fig, (axs) = plt.subplots(nrows=3, ncols=3, figsize=(8,8))
	fig.suptitle(r"$\beta=%lf$" %beta[b])
	axs[0,0].set_xlabel("t")
	axs[0,0].set_ylabel(r"$\beta E(t)/V$")
	axs[0,0].plot(E,'.', markersize=0.05)
	axs[0,2].set_xlabel("t")
	axs[0,2].set_ylabel(r"$A_{\beta E/V}(t)$")
	axs[0,2].plot(acf(E,  nlags=100, fft=True), 'o-') 
	axs[0,1].set_xlabel(r"$\beta E/V$")
	axs[0,1].set_ylabel(r"$P(\beta E/V)$")
	axs[0,1].hist(E, bins=Nbins, density=True)
	axs[0,1].plot(PE.support, PE.density)
	axs[1,0].set_xlabel("t")
	axs[1,0].set_ylabel(r"$U$")
	axs[1,0].plot(U,'.', markersize=0.05)
	axs[1,2].set_xlabel("t")
	axs[1,2].set_ylabel("$A_U(t)$")
	axs[1,2].plot(acf(U,  nlags=100, fft=True), 'o-') 
	axs[1,1].set_xlabel(r"$U$")
	axs[1,1].set_ylabel(r"$P(U)$")
	axs[1,1].hist(U, bins=Nbins, density=True)
	axs[1,1].plot(PU.support, PU.density)
	axs[2,0].set_xlabel("t")
	axs[2,0].set_ylabel("M(t)")
	axs[2,0].plot(M,'.', markersize=0.05)
	axs[2,2].set_xlabel("t")
	axs[2,2].set_ylabel(r"$A_M(t)$")
	axs[2,2].plot(acf(M,  nlags=100, fft=True), 'o-') 
	axs[2,1].set_xlabel("M")
	axs[2,1].set_ylabel("P(M)")
	axs[2,1].hist(M, bins=Nbins, density=True)
	axs[2,1].plot(PM.support, PM.density)

	fig.tight_layout()
	fig.subplots_adjust(top=0.9)
	#plt.show()
	plt.savefig('%s/Probability_Distributions_beta%s.png' %(BASEDIR, beta[b]))
	plt.close()
