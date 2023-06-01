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

def bootstrap_magn(M_data, nblocks, nrs):
    Mboot = np.zeros(nrs)
    Chiboot = np.zeros(nrs)
    Uboot = np.zeros(nrs)
    #bootstrap resampling extract Nblocks with replacement and form a new set of data from which compute Cv, E_err, Cv_err
    for i in range(nrs):
        inds = np.random.randint(nblocks, size=nblocks)
        data = M_data[inds,:]
        data = data.flatten()
        Mboot[i] = np.mean(data)
        Chiboot[i] = np.var(data)
        Uboot[i]= np.mean(np.power(data,4))/(3*np.power(np.mean(np.power(data,2)),2))
    M_err= np.std(Mboot, ddof=1)
    Chi_err= np.std(Chiboot, ddof=1)
    U_err= np.std(Uboot, ddof=1)

    return M_err, Chi_err, U_err, Uboot

nrs=500

U_resampling=np.zeros((nrs, nbeta))

for l in range(len(L)):

    U_mean=np.zeros((nbeta))
    U_err=np.zeros((nbeta))
    M_mean=np.zeros((nbeta))
    M_err=np.zeros((nbeta))

    BASEDIR=("%s/L%d_rho%s_alpha%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s"  %(folder_out, L[l], rho, alpha, eta1, eta2, e,  h, nu, beta_low, beta_high, nMAX, H_init))
    
    data_tau_max=np.loadtxt("%s/tau_max.txt" %BASEDIR, dtype=str)
    tau_max=np.amax(np.array(data_tau_max[1], dtype=float))

    data_transient_time=np.loadtxt("%s/transient_time.txt" %BASEDIR, dtype=str)
    transient_time=int(np.amax(np.array(data_transient_time[1], dtype=float)))

    for b in range(nbeta):
        beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)

        file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
        M=np.asarray(file['Measurements']['m'])

        #cut of the transient regime:
        M=M[transient_time:]
        ####ORDER ACCORDING TO THE RANK#########
        rank=np.asarray(file['Measurements']['rank'])
        rank=rank[transient_time:]
        indices_rank= rank.argsort()
        M=M[indices_rank]
        
        #split the N measurements in Nb blocks according to the autocorrelation time tau
        block_size=int(100*tau_max)
        nblocks=int(len(M)/block_size)

        while((block_size<(100*tau_max)) and (nblocks>20) ):
            nblocks=int(nblocks*0.5)
            block_size=int(len(M)/nblocks)

        M_data=np.zeros((nblocks, block_size))
        for block in range(nblocks):
            M_data[block]=M[block*block_size: (block+1)*block_size]

        M_err, Chi_err, U_err,  U_resampling[:,b]= bootstrap_magn(M_data, nblocks, nrs)

        # meanM_resampling=np.zeros((N_dataset))
        # #bootstrap resampling extract M blocks with replacement and form a new set of data from which compute Cv, E_err, Cv_err
        # for n in range(N_dataset):
        #     resampling= np.random.choice(Nblocks,Nblocks)
        #     M_resampling=M_block[resampling]
        #     meanM_resampling[n]=np.mean(M_resampling)
        #     U_resampling[n, b]=np.mean(np.power(M_resampling,4))/(3*np.power(np.mean(np.power(M_resampling,2)),2))

        # M_mean[b]=np.mean(meanM_resampling)
        M_mean[b]=np.mean(M)
        M_err[b]= M_err
        # M_err[b]=np.std(meanM_resampling)

        # U_mean[b]=np.mean(U_resampling[:,b])
        U_mean[b]=np.mean(np.power(M,4))/(3*np.power(np.mean(np.power(M,2)),2))
        U_err[b]= U_err
        # U_err[b]=np.std(U_resampling[:,b])

    np.savetxt("%s/Binder_alln.txt" %(BASEDIR), U_resampling)

    np.savetxt("%s/Magnetization.txt" %(BASEDIR), (beta, M_mean, M_err))
    np.savetxt("%s/Binder_cumulant.txt" %(BASEDIR), (beta, U_mean, U_err))



