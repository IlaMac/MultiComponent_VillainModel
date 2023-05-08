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
from scipy.optimize import curve_fit

##############################
def linear_fit(x,a,b):
    return a*x +b

##############################

BASEDIR=sys.argv[1]
blow=float(sys.argv[2])
bhigh=float(sys.argv[3])
e=float(sys.argv[4])
h=(sys.argv[5])
nu=float(sys.argv[6])
eta1=(sys.argv[7])
eta2=(sys.argv[8])
rho=sys.argv[9]
alpha=sys.argv[10]
nMAX=sys.argv[11]
Hinit=sys.argv[12]
bmin=float(sys.argv[13])
bmax=float(sys.argv[14])
var=sys.argv[15]
if( (nu).is_integer()): nu=int(nu)
if( (e).is_integer()): e=int(e)

L=[]
for ind in range(16, len(sys.argv)):
    L.append(int(sys.argv[ind]))



####################### U1 CROSSING POINT #######################

HOMEDIR=("%s/L%s_rho%s_alpha%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s" %(BASEDIR, L[0], rho, alpha, eta1, eta2, e, h, nu, blow, bhigh, nMAX, Hinit))
data_Ds=np.loadtxt("%s/Dual_Stiffness.txt" %HOMEDIR)
beta= data_Ds[0]
Ds= data_Ds[1]
err_Ds= data_Ds[2]
start=np.where(beta>bmin)[0][0]
end=np.where(beta<bmax)[0][-1]
nbeta_cut=end-start

N_dataset=100
betac_cross=np.zeros((N_dataset, len(L)-1))

Ds_cross=np.zeros((N_dataset, nbeta_cut, len(L)))
for l in range(len(L)):
    for n in range(N_dataset):
        HOMEDIR=("%s/L%s_rho%s_alpha%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s" %(BASEDIR, L[l], rho, alpha, eta1, eta2, e, h, nu, blow, bhigh, nMAX, Hinit))
        data_Ds=np.loadtxt("%s/Dual_Stiffness.txt" %HOMEDIR)
        beta= data_Ds[0]
        Ds= data_Ds[1]
        err_Ds= data_Ds[2]
        beta=beta[start:end]
        Ds_all=np.loadtxt("%s/DStiffness_alln.txt" %HOMEDIR)
        Ds_cross[n, :, l]=L[l]*Ds_all[n][start:end]

for n in range(N_dataset):   
    for l1 in range(len(L)):
        for l2 in range(l1+1, len(L)):
            #print(beta[np.where(Js_cross[n,:, l1]>Js_cross[n,:, l2])[0][-1]])
            #Return the roots of the (non-linear) equations defined by func(x) = 0 given a starting estimate
            if(len(np.where(Ds_cross[n,:, l1]> (Ds_cross[n,:, l2]))[0])!=0):
                index1=np.where(Ds_cross[n,:, l1]> (Ds_cross[n,:, l2]))[0][0]    
                index2=index1-1
                x1= beta[index1]
                y1= Ds_cross[n,index1, l1]
                x2= beta[index2] 
                y2=Ds_cross[n,index2, l1]
                m1= (y1-y2)/(x1-x2)
                q1= -x2*m1+y2
                y1= Ds_cross[n,index1, l2]
                y2=Ds_cross[n,index2, l2]
                m2= (y1-y2)/(x1-x2)
                q2= -x2*m2+y2
                betac_cross[n,l2-1]=(q2-q1)/(m1-m2)
            else:
                print("nada")

betac_u1_finitesize=[]
err_betac_u1_finitesize=[]

for l2 in range(1, len(L)):
    betac_u1_finitesize.append(np.mean(betac_cross[:,l2-1]))
    err_betac_u1_finitesize.append(np.sqrt(N_dataset-1)*np.std(betac_cross[:,l2-1]))

pair_l=[]

pair_l=np.zeros(len(L)-1)
for l in range(len(L)-1):
    pair_l[l]=1./np.sqrt(L[l]*L[l+1])

np.savetxt("%s/Betac_U1_finitesize_rho%s_eta2%s_e%s_nu%s.txt" %(BASEDIR, rho, eta2, e, nu), (pair_l, np.asarray(betac_u1_finitesize),np.asarray(err_betac_u1_finitesize)))
#print("mean beta_c U(1):", np.mean(betac_u1_finitesize), np.std(betac_u1_finitesize)/np.sqrt(len(betac_u1_finitesize)-1))
plt.errorbar(pair_l, betac_u1_finitesize,yerr= err_betac_u1_finitesize, fmt="o-")

u1_popt, u1_pcov = curve_fit(linear_fit, pair_l , betac_u1_finitesize, sigma=err_betac_u1_finitesize, absolute_sigma=True)
print( var, u1_popt[1], u1_pcov[1,1])

plt.plot(pair_l,linear_fit(pair_l, *u1_popt), ls="--", c="gray" )
plt.show()
