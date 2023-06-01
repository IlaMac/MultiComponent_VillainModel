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

print(sys.argv)


def bootstrap_energy(E_data, nblocks, nrs):
    Eboot = np.zeros(nrs)
    Cboot = np.zeros(nrs)
    #bootstrap resampling extract Nblocks with replacement and form a new set of data from which compute Cv, E_err, Cv_err
    for i in range(nrs):
        inds = np.random.randint(nblocks, size=nblocks)
        data = E_data[inds,:]
        data = data.flatten()
        Eboot[i] = np.mean(data)
        Cboot[i] = np.var(data)
    E_err= np.std(Eboot, ddof=1)
    C_err= np.std(Cboot, ddof=1)

    return E_err, C_err

nrs=500

for l in range(len(L)):

    Cv_mean=np.zeros((nbeta))
    Cv_err=np.zeros((nbeta))
    E_mean=np.zeros((nbeta))
    E_err=np.zeros((nbeta))


    BASEDIR=("%s/L%d_rho%s_alpha%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s"  %(folder_out, L[l], rho, alpha, eta1, eta2, e,  h, nu, beta_low, beta_high, nMAX, H_init))

    data_tau_max=np.loadtxt("%s/tau_max.txt" %BASEDIR, dtype=str)
    tau_max=(np.array(data_tau_max[1][0], dtype=float))

    data_transient_time=np.loadtxt("%s/transient_time.txt" %BASEDIR, dtype=str)
    transient_time=int(np.amax(np.array(data_transient_time[1], dtype=float)))

    for b in range(nbeta):
        beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)
        
        file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
        E=np.asarray(file['Measurements']['U'])

        #cut of the transient regime:
        E=E[transient_time:]
        ####ORDER ACCORDING TO THE RANK#########
        rank=np.asarray(file['Measurements']['rank'])
        rank=rank[transient_time:]
        indices_rank= rank.argsort()
        E=E[indices_rank]
        #split the N measurements in Nblocks blocks according to the autocorrelation time tau

        block_size=int(100*tau_max)
        nblocks=int(len(E)/block_size)
        while((block_size<(100*tau_max)) and (nblocks>20) ):
            Nblocks=int(nblocks*0.5)
            block_size=int(len(E)/nblocks)


        E_data=np.zeros((nblocks, block_size))
        for block in range(nblocks):
            E_data[block]=E[block*block_size: (block+1)*block_size]

        Eerr, Cerr = bootstrap_energy(E_data, nblocks, nrs)

        E_mean[b]=np.mean(E)/(L**3)
        E_err[b]=Eerr/(L**3)
       
        Cv_mean[b]=np.var(E)*(beta[b]*beta[b])/(L**3) 
        Cv_err[b]=Cerr*(beta[b]*beta[b])/(L**3) 


    data_energy=np.vstack((beta, E_mean, E_err))
    np.savetxt("%s/Energy.txt" %(BASEDIR), data_energy)
    np.savetxt("%s/Specific_Heat.txt" %(BASEDIR), (beta, Cv_mean, Cv_err))




