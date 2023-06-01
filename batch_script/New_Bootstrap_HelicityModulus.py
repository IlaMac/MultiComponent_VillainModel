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


NC=2 #number of components
h3=float(h)**3

nrs=500

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

    J1=np.zeros((nrs, nbeta))
    J2=np.zeros((nrs, nbeta))
    Jmixed=np.zeros((nrs, nbeta))
    Jsum=np.zeros((nrs, nbeta))
    Jdiff=np.zeros((nrs, nbeta))

    err_J1=np.zeros((nrs, nbeta))
    err_J2=np.zeros((nrs, nbeta))
    err_Jmixed=np.zeros((nrs, nbeta))
    err_Jsum=np.zeros((nrs, nbeta))
    err_Jdiff=np.zeros((nrs, nbeta))


    BASEDIR=("%s/L%d_rho%s_alpha%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s"  %(folder_out, L[l], rho, alpha, eta1, eta2, e,  h, nu, beta_low, beta_high, nMAX, H_init))

    print(BASEDIR)
    data_tau_max=np.loadtxt("%s/tau_max.txt" %BASEDIR, dtype=str)
    tau_max=np.amax(np.array(data_tau_max[1], dtype=float))
    
    data_transient_time=np.loadtxt("%s/transient_time.txt" %BASEDIR, dtype=str)
    transient_time=int(np.amax(np.array(data_transient_time[1], dtype=float)))

    mean_d1_resampling=np.zeros((nrs, nbeta))
    mean_d1d1_resampling=np.zeros((nrs, nbeta))
    mean_d2_resampling=np.zeros((nrs, nbeta))
    mean_d2d2_resampling=np.zeros((nrs, nbeta))
    mean_d11_resampling=np.zeros((nrs, nbeta))
    mean_d22_resampling=np.zeros((nrs, nbeta))
    mean_d1d2_resampling=np.zeros((nrs, nbeta))
    mean_d12_resampling=np.zeros((nrs, nbeta))

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


        block_size=int(100*tau_max)
        nblocks=int(len(d1)/block_size)
        while((block_size<(100*tau_max)) and (nblocks>20) ):
            nblocks=int(nblocks*0.5)
            block_size=int(len(d1)/nblocks)

        d1_data=np.zeros((nblocks, block_size))
        d2_data=np.zeros((nblocks, block_size))
        d11_data=np.zeros((nblocks, block_size))
        d22_data=np.zeros((nblocks, block_size))
        d12_data=np.zeros((nblocks, block_size))

        for block in range(nblocks):
            d1_data[block]=d1[block*block_size: (block+1)*block_size]
            d2_data[block]=d2[block*block_size: (block+1)*block_size]
            d11_data[block]=d11[block*block_size: (block+1)*block_size]
            d22_data[block]=d22[block*block_size: (block+1)*block_size]
            d12_data[block]=d12[block*block_size: (block+1)*block_size]

        #bootstrap resampling extract Nblocks with replacement and form a new set of data from which compute Cv, E_err, Cv_err
        for i in range(nrs):
            inds = np.random.randint(nblocks, size=nblocks)
            d1_data_i = d1_data[inds,:]
            d1_data_i = d1_data_i.flatten()
            d11_data_i = d11_data[inds,:]
            d11_data_i = d11_data_i.flatten()

            d2_data_i = d2_data[inds,:]
            d2_data_i = d2_data_i.flatten()
            d22_data_i = d22_data[inds,:]
            d22_data_i = d22_data_i.flatten()

            d12_data_i = d12_data[inds,:]
            d12_data_i = d12_data_i.flatten()
            
            mean_d1_resampling[i,b]=np.mean(d1_data_i)
            mean_d1d1_resampling[i,b]=np.mean(d1_data_i*d1_data_i)
            mean_d11_resampling[i,b]=np.mean(d11_data_i)

            mean_d2_resampling[i,b]=np.mean(d2_data_i)
            mean_d2d2_resampling[i,b]=np.mean(d2_data_i*d2_data_i)
            mean_d22_resampling[i,b]=np.mean(d22_data_i)

            mean_d1d2_resampling[i,b]=np.mean(d1_data_i*d2_data_i)
            mean_d12_resampling[i,b]=np.mean(d12_data_i)
          
            # < \Delta H (alpha, beta)/ \Delta \deta_{alpha} >

            J1[i,b]=1./(N) *(mean_d11_resampling[i, b] - (beta[b])*(mean_d1d1_resampling[i, b] - (mean_d1_resampling[i,b]*mean_d1_resampling[i,b]) ) )
            
            J2[i,b]=1./(N) *(mean_d22_resampling[i, b] - (beta[b])*(mean_d2d2_resampling[i, b] - (mean_d2_resampling[i,b]*mean_d2_resampling[i,b]) ) )

            Jmixed[i,b]=1./(N) *(mean_d12_resampling[i, b] - (beta[b])*(mean_d1d2_resampling[i, b] - (mean_d1_resampling[i,b]*mean_d2_resampling[i,b]) ) )
            
            Jsum[i,b]= J1[i,b] +J2[i,b] +2*Jmixed[i,b]

            Jdiff[i,b]= J1[i,b] +J2[i,b] -2*Jmixed[i,b]


        J1_mean[b]=1./(N) *( np.mean(d11) - (beta[b])*(np.mean(d1**2) - np.mean(d1)**2 ) )
        J1_err[b]=(np.std(J1[:,b]))

        J2_mean[b]=1./(N) *( np.mean(d22) - (beta[b])*(np.mean(d2**2) - np.mean(d2)**2 ) )
        J2_err[b]=(np.std(J2[:,b]))

        Jmixed_mean[b]=1./(N) *( np.mean(d12) - (beta[b])*(np.mean(d1*d2) - np.mean(d1)*np.mean(d2) ) )
        Jmixed_err[b]=(np.std(Jmixed[:,b]))

        Jsum_mean[b]= J1_mean[b]+ J2_mean[b] +2*Jmixed_mean[b]
        Jsum_err[b]=(np.std(Jsum[:,b]))

        Jdiff_mean[b]= J1_mean[b]+ J2_mean[b] -2*Jmixed_mean[b]
        Jdiff_err[b]=(np.std(Jdiff[:,b]))


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




