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


folder_out=sys.argv[1]
beta_low=float(sys.argv[2])
beta_high=float(sys.argv[3])
nbeta=int(sys.argv[4])
e=float(sys.argv[5])
h=(sys.argv[6])
nu=float(sys.argv[7])
eta1=(sys.argv[8])
eta2=(sys.argv[9])
rho=sys.argv[10]
alpha=sys.argv[11]
nMAX=sys.argv[12]
H_init=sys.argv[13]
flag_init=(sys.argv[14])

if( (nu).is_integer()): nu=int(nu)
if( (e).is_integer()): e=int(e)
#if( (h).is_integer()): h=int(h)


L=[]
for ind in range(15, len(sys.argv)):
    L.append(int(sys.argv[ind]))

beta=np.zeros((nbeta))

color=iter(plt.cm.rainbow(np.linspace(0,1,len(L)+1)))


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{bm}')
fig, ax1 = plt.subplots(nrows=2, ncols=2, figsize=(4,8))
fig.suptitle(r"$h=%s$; $e=%s$; $\nu=%s$" %(h, e, nu))
ax1[0,0].set_xlabel(r"$\beta$")
ax1[0,0].set_ylabel(r"$M_1$")
ax1[1,0].set_xlabel(r"$\beta$")
ax1[1,0].set_ylabel(r"$U_2(\phi_1)$")
ax1[0,1].set_xlabel(r"$\beta$")
ax1[0,1].set_ylabel(r"$M_2$")
ax1[1,1].set_xlabel(r"$\beta$")
ax1[1,1].set_ylabel(r"$U_2(\phi_2)$")

for l in range(len(L)):
    U_1_mean=np.zeros((nbeta))
    U_1_err=np.zeros((nbeta))
    M_1_mean=np.zeros((nbeta))
    M_1_err=np.zeros((nbeta))
    
    U_2_mean=np.zeros((nbeta))
    U_2_err=np.zeros((nbeta))
    M_2_mean=np.zeros((nbeta))
    M_2_err=np.zeros((nbeta))
    
    N_dataset=100

    meanM1_resampling=np.zeros((N_dataset, nbeta)) 
    U1_resampling=np.zeros((N_dataset, nbeta))     
    
    meanM2_resampling=np.zeros((N_dataset, nbeta))
    U2_resampling=np.zeros((N_dataset, nbeta))
    
    BASEDIR=("%s/L%d_rho%s_alpha%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s"  %(folder_out, L[l], rho, alpha, eta1, eta2, e,  h, nu, beta_low, beta_high, nMAX, H_init))

    c_m=next(color)
    
    data_tau_max=np.loadtxt("%s/tau_max.txt" %BASEDIR, dtype=str)
    tau_max=np.amax(np.array(data_tau_max[1], dtype=float))

    data_transient_time=np.loadtxt("%s/transient_time.txt" %BASEDIR, dtype=str)
    transient_time=int(np.amax(np.array(data_transient_time[1], dtype=float)))


    for b in range(nbeta):
        beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)

        file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
        M=np.asarray(file['Measurements']['m_phase'])

        #cut of the transient regime:
      	#In my c++ code I have computed M^2= Mx^2 + My^2
        M1=np.sqrt(M[transient_time:,0])
        M2=np.sqrt(M[transient_time:,1])
        #split the N measurements in Nb blocks according to the autocorrelation time tau
        Nblocks=100
        block_size=int(len(M1)/Nblocks)
        while((block_size<(20*tau_max)) and (Nblocks>20) ):
            Nblocks=int(Nblocks*0.5)
            block_size=int(len(M1)/Nblocks)
        M1_block=np.zeros((Nblocks, block_size))
        M2_block=np.zeros((Nblocks, block_size))

        for block in range(Nblocks):
            M1_block[block]=M1[block*block_size: (block+1)*block_size]
            M2_block[block]=M2[block*block_size: (block+1)*block_size]

        #bootstrap resampling extract M blocks with replacement and form a new set of data from which compute Cv, E_err, Cv_err
        for n in range(N_dataset):
            resampling= np.random.choice(Nblocks,Nblocks)
            M1_resampling=M1_block[resampling]
            meanM1_resampling[n,b]=np.mean(M1_resampling)
            #definition for 2 component (XY like) order parameter
            U1_resampling[n,b]=2*(1- np.mean(np.power(M1_resampling,4))/(2*np.power(np.mean(np.power(M1_resampling,2)),2)))
            
            M2_resampling=M2_block[resampling]
            meanM2_resampling[n,b]=np.mean(M2_resampling)
            #definition for 2 component (XY like) order parameter
            U2_resampling[n,b]=2*(1- np.mean(np.power(M2_resampling,4))/(2*np.power(np.mean(np.power(M2_resampling,2)),2)))

        M_1_mean[b]=np.mean(meanM1_resampling[:, b])
        M_1_err[b]=np.std(meanM1_resampling[:, b])
        U_1_mean[b]=np.mean(U1_resampling[:, b])
        U_1_err[b]=np.std(U1_resampling[:, b])

        M_2_mean[b]=np.mean(meanM1_resampling[:, b])
        M_2_err[b]=np.std(meanM1_resampling[:, b])
        U_2_mean[b]=np.mean(U1_resampling[:, b])
        U_2_err[b]=np.std(U1_resampling[:, b])


    np.savetxt("%s/Magnetization_phase.txt" %(BASEDIR), (beta, M_1_mean, M_1_err,  M_1_mean, M_2_err))
    np.savetxt("%s/Binder_phase_cumulant.txt" %(BASEDIR), (beta, U_1_mean, U_1_err, U_2_mean, U_2_err))

    ax1[0,0].plot(beta, M_1_mean, '-', c=c_m)
    ax1[0,0].errorbar(beta, M_1_mean, yerr=M_1_err, capsize=2, c=c_m,label="L=%s" %L[l])
    ax1[1,0].errorbar(beta, U_1_mean, yerr=U_1_err, c=c_m, capsize=2)
    ax1[0,1].plot(beta, M_2_mean, '-', c=c_m)
    ax1[0,1].errorbar(beta, M_2_mean, yerr=M_2_err, capsize=2, c=c_m,label="L=%s" %L[l])
    ax1[1,1].errorbar(beta, U_2_mean, yerr=U_2_err, c=c_m, capsize=2)

ax1[0,0].legend(loc="best")
fig.tight_layout()
fig.savefig("%s/Magnetization_phase_h%s_bmin%s_bmax%s.png" %(folder_out, h, beta_low, beta_high))
#plt.show()


