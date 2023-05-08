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
import h5py

beta_low=float(sys.argv[1])
beta_high=float(sys.argv[2])
nbeta=int(sys.argv[3])
BASEDIR=sys.argv[4]
L=sys.argv[5]
nu=float(sys.argv[6])
e=float(sys.argv[7])
eta= float(sys.argv[8])
beta=np.zeros((nbeta))


   
Observables=np.array(["E", "m", "m_phase", "D2H_Dd2i", "DH_Ddi"])    

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{bm}')

tau_max=np.zeros((len(Observables)))

#Observables=np.array(["D2H_Dd2i", "DH_Ddi"])
for name in range(len(Observables)):
    print(L, Observables[name]
    )

    tau=np.zeros((nbeta))

    for b in range(nbeta):
        beta[b]=beta_low +b*(beta_high -beta_low)/(nbeta-1)
        Obs_mean=np.zeros((nbeta))
        Obs_var=np.zeros((nbeta))
        file=h5py.File('%s/beta_%d/Output.h5' %(BASEDIR, b), 'r')
        Obs=np.asarray(file['Measurements']['%s' %(Observables[name])])

        ######################################
        if((Observables[name]=="D2H_Dd2i") or (Observables[name]=="DH_Ddi") or (Observables[name]=="m_phase")):
            #print(np.shape(Obs), Obs)
            Obs=Obs[:,0]

        ####ORDER ACCORDING TO THE RANK#########
        rank=np.asarray(file['Measurements']['rank'])
        indices_rank= rank.argsort()
        old_Obs=Obs
        Obs=Obs[indices_rank]

        A_Obs=acf(Obs, nlags=int(len(Obs)/10), fft=True)

        #plt.plot(A_Obs)
        #plt.show()
        temp_tau=[]
        time_int=10
        while(time_int<= len(A_Obs)):
            temp_tau=np.append(temp_tau, np.sum(A_Obs[:time_int]))
            time_int=time_int+100
        tau[b]=np.amax(temp_tau)


    tau_max[name]=np.amax(tau)

data=np.vstack((Observables, tau_max))

print(data)
np.savetxt("%s/tau_max.txt" %BASEDIR, data, fmt="%s")




