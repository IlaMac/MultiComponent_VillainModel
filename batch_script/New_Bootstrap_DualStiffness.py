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
nMAX=sys.argv[11]
H_init=sys.argv[12]
flag_init=(sys.argv[13])

if( (nu).is_integer()): nu=int(nu)
if( (e).is_integer()): e=int(e)
#if( (h).is_integer()): h=int(h)


L=[]
for ind in range(14, len(sys.argv)):
    L.append(int(sys.argv[ind]))

beta=np.zeros((nbeta))

print(sys.argv)


######################################################################################################
#  ! My definition of rho_z differs from the Daniel's one by a factor h^2 which I add here. #
######################################################################################################

factor=float(h)**2

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=18)
plt.rc('text.latex', preamble=r'\usepackage{bm}')
fig, ax1 = plt.subplots(1, 1, figsize=(6,6))
ax1.set_title(r"$h=%s$; $e=%s$; $\nu=%s$" %(h, e, nu))
ax1.set_xlabel(r"$\beta$")
ax1.set_ylabel(r"$L\rho$")

color=iter(plt.cm.rainbow(np.linspace(0,1,len(L)+1)))


for l in range(len(L)):

    Ds_mean=np.zeros((nbeta))
    Ds_err=np.zeros((nbeta))

    c_m=next(color)

    N_dataset=100

    Ds_all=np.zeros((N_dataset, nbeta))

    BASEDIR=("%s/L%d_rho%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s"  %(folder_out, L[l], rho, eta1, eta2, e,  h, nu, beta_low, beta_high, nMAX, H_init))

    data_tau_max=np.loadtxt("%s/tau_max.txt" %BASEDIR, dtype=str)
    tau_max=np.amax(np.array(data_tau_max[1], dtype=float))

    data_transient_time=np.loadtxt("%s/transient_time.txt" %BASEDIR, dtype=str)
    transient_time=int(np.amax(np.array(data_transient_time[1], dtype=float)))

    for b in range(nbeta):
        beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)

        file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
        Ds=np.asarray(file['Measurements']['ds'])

        #cut of the transient regime:
        Ds=factor*Ds[transient_time:]
        ####ORDER ACCORDING TO THE RANK#########
        rank=np.asarray(file['Measurements']['rank'])
        rank=rank[transient_time:]
        indices_rank= rank.argsort()
        Ds=Ds[indices_rank]
        
        #split the N measurements in Nblocks blocks according to the autocorrelation time tau
        Nblocks=100
        block_size=int(len(Ds)/Nblocks)
        while((block_size<(20*tau_max)) and (Nblocks>20) ):
            Nblocks=int(Nblocks*0.5)
            block_size=int(len(Ds)/Nblocks)

        Ds_block=np.zeros((Nblocks, block_size))

        for block in range(Nblocks):
            Ds_block[block]=Ds[block*block_size: (block+1)*block_size]

        meanDs_resampling=np.zeros((N_dataset))
        #bootstrap resampling extract M blocks with replacement and form a new set of data from which compute Cv, E_err, Cv_err
        for n in range(N_dataset):
            resampling= np.random.choice(Nblocks,Nblocks)
            Ds_resampling=Ds_block[resampling]
            meanDs_resampling[n]=np.mean(Ds_resampling)

            Ds_all[n,b]=np.mean(Ds_resampling)

        Ds_mean[b]=np.mean(meanDs_resampling)
        Ds_err[b]=np.sqrt(N_dataset-1)*np.std(meanDs_resampling)

    np.savetxt("%s/DStiffness_alln.txt" %(BASEDIR), Ds_all)
    np.savetxt("%s/Dual_Stiffness.txt" %(BASEDIR), (beta, Ds_mean, Ds_err))


#     ax1.plot(beta, Ds_mean, '-', c=c_m)
#     ax1.errorbar(beta, Ds_mean, yerr=Ds_err, capsize=2, c=c_m, label="L=%s" %L[l])

# ax1.legend(loc="best")
# plt.tight_layout()
# plt.savefig("%s/Dual_Stiffness_h%s_bmin%s_bmax%s.png" %(folder_out, h, beta_low, beta_high))

