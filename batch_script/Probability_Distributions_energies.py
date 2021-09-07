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


HOMEDIR=sys.argv[1]
blow=float(sys.argv[2])
bhigh=float(sys.argv[3])
e=float(sys.argv[4])
h=int(sys.argv[5])
nu=float(sys.argv[6])
eta1=(sys.argv[7])
eta2=(sys.argv[8])
rho=sys.argv[9]
alpha=sys.argv[10]
nMAX=sys.argv[11]
Hinit=sys.argv[12]
nbetas=int(sys.argv[13])
if( (nu).is_integer()): nu=int(nu)
if( (e).is_integer()): e=int(e)
LLIST=[]
for ind in range(14, len(sys.argv)):
    LLIST.append(int(sys.argv[ind]))


beta=np.zeros((nbetas))

Nbins=50
for L in LLIST:
	L=int(L)
	BASEDIR=("%s/L%s_rho%s_alpha%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s" %(HOMEDIR, L, rho, alpha, eta1, eta2, e, h, nu, blow, bhigh, nMAX, Hinit))
	folder_distributions=("%s/PDistributions" %HOMEDIR)
	try:
	    os.makedirs(folder_distributions)    
	    print("Directory " , folder_distributions ,  " Created ")
	except FileExistsError:
	    print("Directory " , folder_distributions ,  " already exists") 
	temp_length=[]

	folder_distributions_L=("%s/PDistributions/L_%s_bmin%s_bmax%s" %(HOMEDIR,L, blow, bhigh))
	try:
	    os.makedirs(folder_distributions_L)    
	    print("Directory " , folder_distributions_L ,  " Created ")
	except FileExistsError:
	    print("Directory " , folder_distributions_L ,  " already exists") 
	temp_length=[]
	for b in range(nbetas):
		beta[b]=blow +b*(bhigh -blow)/(nbetas-1)

		file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')

		E=np.asarray(file['Measurements']['U'])
		E=E/pow((L*h),3)
		histE, bin_edgesE = np.histogram(E)
		PE = sm.nonparametric.KDEUnivariate(E)
		PE.fit()

		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
		plt.rc('text.latex', preamble=r'\usepackage{bm}')
		fig, (axs) = plt.subplots(nrows=3, ncols=1, figsize=(6,6))
		fig.suptitle(r"$\beta=%lf$" %beta[b])
		axs[0].set_xlabel("t")
		axs[0].set_ylabel(r"$E_{TOT}(t)/V$")
		axs[0].plot(E,'.', markersize=0.05)
		axs[2].set_xlabel("t")
		axs[2].set_ylabel(r"$A_{E_{TOT}/V}(t)$")
		axs[2].plot(acf(E,  nlags=100, fft=True), 'o-') 
		axs[1].set_xlabel(r"$E_{TOT}/V$")
		axs[1].set_ylabel(r"$P( E_{TOT}/V)$")
		axs[1].hist(E, bins=Nbins, density=True)
		axs[1].plot(PE.support, PE.density)
		fig.tight_layout()
		fig.subplots_adjust(top=0.9)
		#plt.show()
		plt.savefig('%s/Probability_Distributions_Energy_beta%s.png' %(folder_distributions_L, b))
		plt.close()
