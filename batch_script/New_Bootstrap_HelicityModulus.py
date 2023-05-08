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


NC=2 #number of components
h3=float(h)**3


for l in range(len(L)):
    
    N=L[l]**3

    Jsum_mean=np.zeros((nbeta))
    Jsum_err=np.zeros((nbeta))

    Jdiff_mean=np.zeros((nbeta))
    Jdiff_err=np.zeros((nbeta))

    J1_mean=np.zeros((nbeta))
    J1_err=np.zeros((nbeta))

    J2_mean=np.zeros((nbeta))
    J2_err=np.zeros((nbeta))

    Jmixed_mean=np.zeros((nbeta))
    Jmixed_err=np.zeros((nbeta))

    N_dataset=100

    J1=np.zeros((N_dataset, nbeta))
    J2=np.zeros((N_dataset, nbeta))
    Jmixed=np.zeros((N_dataset, nbeta))
    Jsum=np.zeros((N_dataset, nbeta))
    Jdiff=np.zeros((N_dataset, nbeta))

    err_J1=np.zeros((N_dataset, nbeta))
    err_J2=np.zeros((N_dataset, nbeta))
    err_Jmixed=np.zeros((N_dataset, nbeta))
    err_Jsum=np.zeros((N_dataset, nbeta))
    err_Jdiff=np.zeros((N_dataset, nbeta))


    BASEDIR=("%s/L%d_rho%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s"  %(folder_out, L[l], rho, eta1, eta2, e,  h, nu, beta_low, beta_high, nMAX, H_init))

    print(BASEDIR)
    data_tau_max=np.loadtxt("%s/tau_max.txt" %BASEDIR, dtype=str)
    tau_max=np.amax(np.array(data_tau_max[1], dtype=float))
    
    data_transient_time=np.loadtxt("%s/transient_time.txt" %BASEDIR, dtype=str)
    transient_time=int(np.amax(np.array(data_transient_time[1], dtype=float)))

    mean_d1_resampling=np.zeros((N_dataset, nbeta))
    mean_d1d1_resampling=np.zeros((N_dataset, nbeta))
    mean_d2_resampling=np.zeros((N_dataset, nbeta))
    mean_d2d2_resampling=np.zeros((N_dataset, nbeta))
    mean_d11_resampling=np.zeros((N_dataset, nbeta))
    mean_d22_resampling=np.zeros((N_dataset, nbeta))
    mean_d1d2_resampling=np.zeros((N_dataset, nbeta))
    mean_d12_resampling=np.zeros((N_dataset, nbeta))

    for b in range(nbeta):
        beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)

        file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')

        DH_Ddi=np.asarray(file['Measurements']['DH_Ddi'])
        D2H_Dd2i=np.asarray(file['Measurements']['D2H_Dd2i'])
        D2H_Dd12=np.asarray(file['Measurements']['D2H_Dd12'])


        d1=DH_Ddi[:,0]
        d2=DH_Ddi[:,1]
        d11=D2H_Dd2i[:,0]
        d22=D2H_Dd2i[:,1]
        d12=D2H_Dd12

        d1=d1[transient_time:]
        d2=d2[transient_time:]
        d11=d11[transient_time:]
        d22=d22[transient_time:]
        d12=d12[transient_time:]

        ####ORDER ACCORDING TO THE RANK#########
        rank=np.asarray(file['Measurements']['rank'])
        rank=rank[transient_time:]
        indices_rank= rank.argsort()
        d1=d1[indices_rank]
        d2=d2[indices_rank]
        d11=d11[indices_rank]
        d22=d22[indices_rank]
        d12=d12[indices_rank]

        Nblocks=100
        block_size=int(len(d1)/Nblocks)
        while((block_size<(20*tau_max)) and (Nblocks>20) ):
            Nblocks=int(Nblocks*0.5)
            block_size=int(len(d1)/Nblocks)

        d1_block=np.zeros((Nblocks, block_size))
        d2_block=np.zeros((Nblocks, block_size))
        d11_block=np.zeros((Nblocks, block_size))
        d22_block=np.zeros((Nblocks, block_size))
        d12_block=np.zeros((Nblocks, block_size))

        for block in range(Nblocks):
            d1_block[block]=d1[block*block_size: (block+1)*block_size]
            d2_block[block]=d2[block*block_size: (block+1)*block_size]
            d11_block[block]=d11[block*block_size: (block+1)*block_size]
            d22_block[block]=d22[block*block_size: (block+1)*block_size]
            d12_block[block]=d12[block*block_size: (block+1)*block_size]

        #bootstrap resampling extract Nblocks with replacement and form a new set of data from which compute Cv, E_err, Cv_err
        for n in range(N_dataset):
            resampling= np.random.choice(Nblocks,Nblocks)
            d1_resampling=d1_block[resampling]
            d1d1_resampling=d1_block[resampling]*d1_block[resampling]
            d11_resampling=d11_block[resampling]
            
            d2_resampling=d2_block[resampling]
            d2d2_resampling=d2_block[resampling]*d2_block[resampling]
            d22_resampling=d22_block[resampling]

            d1d2_resampling=d1_block[resampling]*d2_block[resampling]
            d12_resampling=d12_block[resampling]
            
            mean_d1_resampling[n,b]=np.mean(d1_resampling)
            mean_d1d1_resampling[n,b]=np.mean(d1d1_resampling)
            mean_d11_resampling[n,b]=np.mean(d11_resampling)

            mean_d2_resampling[n,b]=np.mean(d2_resampling)
            mean_d2d2_resampling[n,b]=np.mean(d2d2_resampling)
            mean_d22_resampling[n,b]=np.mean(d22_resampling)
            

            mean_d1d2_resampling[n,b]=np.mean(d1d2_resampling)
            mean_d12_resampling[n,b]=np.mean(d12_resampling)
          
            # < \Delta H (alpha, beta)/ \Delta \deta_{alpha} >

            J1[n,b]=1./(N) *(mean_d11_resampling[n, b] - (beta[b])*(mean_d1d1_resampling[n, b] - (mean_d1_resampling[n,b]*mean_d1_resampling[n,b]) ) )
            
            J2[n,b]=1./(N) *(mean_d22_resampling[n, b] - (beta[b])*(mean_d2d2_resampling[n, b] - (mean_d2_resampling[n,b]*mean_d2_resampling[n,b]) ) )

            Jmixed[n,b]=1./(N) *(mean_d12_resampling[n, b] - (beta[b])*(mean_d1d2_resampling[n, b] - (mean_d1_resampling[n,b]*mean_d2_resampling[n,b]) ) )
            
            Jsum[n,b]= J1[n,b] +J2[n,b] +2*Jmixed[n,b]

            Jdiff[n,b]= J1[n,b] +J2[n,b] -2*Jmixed[n,b]


        J1_mean[b]=np.mean(J1[:,b])
        J1_err[b]=(np.sqrt(N_dataset-1)*np.std(J1[:,b]))

        J2_mean[b]=np.mean(J2[:,b])
        J2_err[b]=(np.sqrt(N_dataset-1)*np.std(J2[:,b]))

        Jmixed_mean[b]=np.mean(Jmixed[:,b])
        Jmixed_err[b]=(np.sqrt(N_dataset-1)*np.std(Jmixed[:,b]))

        Jsum_mean[b]=np.mean(Jsum[:,b])
        Jsum_err[b]=(np.sqrt(N_dataset-1)*np.std(Jsum[:,b]))

        Jdiff_mean[b]=np.mean(Jdiff[:,b])
        Jdiff_err[b]=(np.sqrt(N_dataset-1)*np.std(Jdiff[:,b]))



    np.savetxt("%s/Jsum_alln.txt" %(BASEDIR), Jsum)
    np.savetxt("%s/Jdiff_alln.txt" %(BASEDIR), Jdiff)


    data_helicity=np.array([beta, Jsum_mean, Jsum_err])
    data_helicity=np.transpose(data_helicity)
    np.savetxt("%s/Helicity_modulus_sum.txt" %(BASEDIR), data_helicity, fmt=['%lf','%lf', '%lf'])


    data_helicity=np.array([beta, Jdiff_mean, Jdiff_err])
    data_helicity=np.transpose(data_helicity)
    np.savetxt("%s/Helicity_modulus_diff.txt" %(BASEDIR), data_helicity, fmt=['%lf','%lf', '%lf'])

    data_helicity_singlec=np.array([beta, J1_mean, J1_err, J2_mean, J2_err])
    data_helicity_singlec=np.transpose(data_helicity_singlec)
    np.savetxt("%s/Helicity_modulus_singlecomponents.txt" %(BASEDIR), data_helicity_singlec, fmt=['%lf','%lf', '%lf','%lf', '%lf'])

    data_helicity_mixedc=np.array([beta, Jmixed_mean, Jmixed_err])
    data_helicity_mixedc=np.transpose(data_helicity_mixedc)
    np.savetxt("%s/Helicity_modulus_mixedcomponents.txt" %(BASEDIR), data_helicity_mixedc, fmt=['%lf','%lf', '%lf'])



# color=iter(plt.cm.rainbow(np.linspace(0,1,len(L)+1)))
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# plt.rc('text.latex', preamble=r'\usepackage{bm}')
# fig, ax1 = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(8,8))
# fig.suptitle(r"$h=%s$; $e=%s$; $\nu=%s$" %(h, e, nu))
# ax1[0,0].set_xlabel(r"$\beta$")
# ax1[0,0].set_ylabel(r"$\Upsilon_1$")
# ax1[0,1].set_xlabel(r"$\beta$")
# ax1[0,1].set_ylabel(r"$\Upsilon_2$")
# ax1[1,0].set_xlabel(r"$\beta$")
# ax1[1,0].set_ylabel(r"$\Upsilon_{sum}$")
# ax1[1,1].set_xlabel(r"$\beta$")
# ax1[1,1].set_ylabel(r"$\Upsilon_{diff}$")
# color=iter(plt.cm.rainbow(np.linspace(0,1,len(L)+1)))
# c_m=next(color)
for l in range(len(L)):
    # c_m=next(color)

    data_J= np.loadtxt("%s/Helicity_modulus_singlecomponents.txt" %BASEDIR)
    beta= data_J[:, 0]
    J1= data_J[:,1]
    err_J1= data_J[:, 2]
    J2= data_J[:, 3]
    err_J2= data_J[:, 4]  
    beta, J_sum, err_J_sum=np.loadtxt("%s/Helicity_modulus_sum.txt" %BASEDIR,  usecols=(0,1,2), unpack=True )
    beta, J_diff, err_J_diff=np.loadtxt("%s/Helicity_modulus_diff.txt" %BASEDIR,  usecols=(0,1,2), unpack=True )
    


#     ax1[0,0].errorbar(beta, L[l]*J1, yerr=err_J1, fmt='o-', c=c_m, label="%s" %L[l])
#     ax1[0,1].errorbar(beta, L[l]*J2, yerr=err_J2, fmt='o-', c=c_m)
#     ax1[1,0].errorbar(beta, J_sum, yerr=err_J_sum, fmt='o-', c=c_m)
#     #ax1[1,1].errorbar(beta, J_diff, yerr=err_J_diff, fmt='o-', c=c_m)
#     for n in range(N_dataset):
#         ax1[1,1].plot(beta, Jdiff[n,:], '-')


# ax1[0,0].legend(loc="best")
# fig.tight_layout()
# fig.savefig("%s/Helicity_h%s_bmin%s_bmax%s.png" %(BASEDIR, h, beta_low, beta_high))
#plt.show()

