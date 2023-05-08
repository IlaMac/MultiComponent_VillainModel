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
import scipy.stats as stats

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.serif'] = 'Computer Modern'
plt.rcParams['axes.linewidth']  = 3.0
plt.rcParams['axes.labelsize']  = 20
plt.rcParams.update({'font.size': 22})
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22
plt.rcParams['xtick.major.size'] = 4
plt.rcParams['ytick.major.size'] = 4
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['legend.fontsize']  = 22
plt.rcParams['legend.frameon']  = True
plt.rcParams["legend.fancybox"] = True
plt.rcParams["legend.shadow"] = True
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["axes.facecolor"] = 'white'
plt.rcParams["axes.edgecolor"] = 'black'
##############################
def linear_fit(x,a,b):
    return a*x +b
##############################
##############################
def exp_fit(x,a,b):
    return b*np.exp(x*a)

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

# blow=0.55
# bhigh=0.65
# e=2.5
# h=1
# nu=0
# eta1=0
# eta2=0.1
# rho=1
# alpha=1
# nMAX=30
# Hinit=1
# bmin=0.6
# bmax=0.65
# var=e


# L=[8, 10, 12, 16, 20, 24, 32]
# BASEDIR=("/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_%s/nu_%s/eta2_%s/h_%s/" %(e, nu, eta2, h))

# ####################### Z2 CROSSING POINT #######################


HOMEDIR=("%s/L%s_rho%s_alpha%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s" %(BASEDIR, L[0], rho, alpha, eta1, eta2, e, h, nu, blow, bhigh, nMAX, Hinit))
# beta, J_diff, err_J_diff=np.loadtxt("%s/Helicity_modulus_diff.txt" %HOMEDIR,  usecols=(0,1,2), unpack=True )
# Js_all=np.loadtxt("%s/Jdiff_alln.txt" %HOMEDIR)
data_U=np.loadtxt("%s/Binder_cumulant.txt" %HOMEDIR)
beta= data_U[0]
U= data_U[1]
err_U= data_U[2]
U_all=np.loadtxt("%s/Binder_alln.txt" %HOMEDIR)
start=np.where(beta>bmin)[0][0]
end=np.where(beta<bmax)[0][-1]
nbeta_cut=end-start

N_dataset=100
betac_cross=np.zeros((N_dataset, len(L)-1))


# Js_cross=np.zeros((N_dataset, nbeta_cut, len(L)))
U_cross=np.zeros((N_dataset, nbeta_cut, len(L)))

for n in range(N_dataset):
    for l in range(len(L)):
        HOMEDIR=("%s/L%s_rho%s_alpha%s_eta1%s_eta2%s_e%s_h%s_nu%s_bmin%s_bmax%s_nMAX%s_init%s" %(BASEDIR, L[l], rho, alpha, eta1, eta2, e, h, nu, blow, bhigh, nMAX, Hinit))
        # beta, J_diff, err_J_diff=np.loadtxt("%s/Helicity_modulus_diff.txt" %HOMEDIR,  usecols=(0,1,2), unpack=True )
        # beta=beta[start:end]
        # Js_all=np.loadtxt("%s/Jdiff_alln.txt" %HOMEDIR)
        # Js_cross[n, :, l]=L[l]*Js_all[n][start:end]

        data_U=np.loadtxt("%s/Binder_cumulant.txt" %HOMEDIR)
        beta= data_U[0]
        U= data_U[1]
        err_U= data_U[2]
        beta=beta[start:end]
        U_all=np.loadtxt("%s/Binder_alln.txt" %HOMEDIR)
        U_cross[n, :, l]=U_all[n][start:end]

betac_z2_finitesize=[]
err_betac_z2_finitesize=[]
for l1 in range(len(L)-1):
    l2=l1+1
    betac_cross_l1=[]
    for n in range(N_dataset):     
        # if(len(np.where(U_cross[n,:, l1]< (U_cross[n,:, l2]))[0])!=0):
        #     index1=np.where(U_cross[n,:, l1]< (U_cross[n,:, l2]))[0][-1]    
        #     index2=index1-1
        condition_1=np.where(U_cross[n,:, l2] <1)[0]
        condition_2=np.where(U_cross[n,:, l2]>(U_cross[n,:, l1]))[0]
        condition=np.intersect1d(condition_1, condition_2)
        #print(condition[-1])
        if(condition.size>0):
            # print(condition)
            index1=condition[-1]
            index2=index1+1
            x1= beta[index1]
            y1= U_cross[n,index1, l1]
            x2= beta[index2] 
            y2=U_cross[n,index2, l1]
            m1= (y1-y2)/(x1-x2)
            q1= -x2*m1+y2
            y1= U_cross[n,index1, l2]
            y2=U_cross[n,index2, l2]
            m2= (y1-y2)/(x1-x2)
            q2= -x2*m2+y2
            betac_test=(q2-q1)/(m1-m2)
            #betac_cross[n,l2-1]=(q2-q1)/(m1-m2)
            if((betac_test>bmin) and (betac_test< bmax)):
                betac_cross_l1.append(betac_test)
        #else:
            #print("nada")

    betac_z2_finitesize.append(np.mean(betac_cross_l1))
    err_betac_z2_finitesize.append(np.sqrt(len(betac_cross_l1)-1)*np.std(betac_cross_l1))

# betac_z2_finitesize=[]
# err_betac_z2_finitesize=[]

# for l2 in range(1, len(L)):
#     betac_z2_finitesize.append(np.mean(betac_cross[:,l2-1]))
#     err_betac_z2_finitesize.append(np.sqrt(N_dataset-1)*np.std(betac_cross[:,l2-1]))

pair_l=[]

pair_l=np.zeros(len(L)-1)
for l in range(len(L)-1):
    pair_l[l]=1./np.sqrt(L[l]*L[l+1])
    
np.savetxt("%s/Betac_z2_finitesize_rho%s_eta2%s_e%s_nu%s.txt" %(BASEDIR, rho, eta2, e, nu), (pair_l, np.asarray(betac_z2_finitesize),np.asarray(err_betac_z2_finitesize)))
plt.rc('text',usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{bm}')
fig, ax1 = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(6,6))
z2_popt, z2_pcov = curve_fit(linear_fit, pair_l , betac_z2_finitesize, sigma=err_betac_z2_finitesize, absolute_sigma=True)

########## GOODNES OF THE FIT #############
residuals = (betac_z2_finitesize)- linear_fit(pair_l, *z2_popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((betac_z2_finitesize-np.mean(betac_z2_finitesize))**2)
r_squared = 1 - (ss_res/ ss_tot)
# print("Z2", ss_res, ss_tot, r_squared)
#perform Chi-Square Goodness of Fit Test
# print(stats.chisquare(f_obs= betac_z2_finitesize, f_exp=linear_fit(pair_l, *z2_popt)))

##############################################
threshold=0.5

if(ss_res < threshold):        
    print( var, z2_popt[1], z2_pcov[1,1])
else:
    err_betac_z2_finitesize=np.asarray(err_betac_z2_finitesize)
    print( var, np.mean(betac_z2_finitesize), np.sqrt(np.sum(err_betac_z2_finitesize*err_betac_z2_finitesize))/len(err_betac_z2_finitesize))

bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)

ax1.errorbar(pair_l, betac_z2_finitesize,yerr= err_betac_z2_finitesize, fmt="o-")
ax1.plot(pair_l,linear_fit(pair_l, *z2_popt), ls="--", c="gray" , label= r"$$RSS=%.2lf$$ $$RSS_{threshold}=%s$$" %(ss_res, threshold))
ax1.set_ylabel(r"$\beta_c(Z_2)$")
ax1.set_xlabel(r"$(L_i L_{i+1})^{-1/2}$")
ax1.text(0.58, 0.1, r"$$\beta_c^{fit}=%.3lf$$ $$\beta_c^{mean}=%.3lf$$ " %(z2_popt[1], np.mean(betac_z2_finitesize)),  transform=ax1.transAxes,  bbox=bbox_props) 
ax1.legend(loc="best")
fig.suptitle(r"$\nu=%s$; $e=%s$; $\eta_2=%s$" %(nu, e, eta2))
fig.tight_layout()
fig.subplots_adjust(top=0.9)
fig.savefig("%s/betacZ2_nu%s_e%s_eta2%s_alpha%s.png" %(BASEDIR, nu, e, eta2, alpha))
#plt.show()
plt.close()
