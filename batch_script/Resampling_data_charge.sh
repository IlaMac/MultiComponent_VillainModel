#!/bin/bash

################# This script computes for a set of L (the array is defined in the .py files): ####################################
#	-the maximum transient time among all the observables
#	-the maximum autocorrelation time among all the observables
#	-Given these two times it perform a bootstrap resampling to compute the mean value and the variance of the observables 
###################################################################################################################################



############# Parameters of the Hamiltonian ##################
H_rho=1
H_alpha=1
H_eta1=0
H_eta2=100
H_e=0
H_h=1
H_nu=0.3
H_blow=0.31
H_bhigh=0.325
nMAX=30
H_init="1"
flag=1 #vecchia nomenclatura cartella senza _init${H_init}

nbeta=64

#LList="\"[[8] [10]]\""

#LList=("8 10 12")
LList=("8 10 12 16 20 24 32")


BASEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_${H_e}/nu_${H_nu}/eta2_${H_eta2}/h_${H_h}"

for L in $LList; do
    DIRECTORY=$BASEDIR/L${L}_rho${H_rho}_alpha${H_alpha}_eta1${H_eta1}_eta2${H_eta2}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_nMAX${nMAX}_init${H_init}
    python3 New_Autocorr_time.py ${H_blow} ${H_bhigh} ${nbeta} ${DIRECTORY} ${L} ${H_nu} ${H_e} ${H_eta1}
    python3 New_LogBoxing.py ${H_blow} ${H_bhigh} ${nbeta} ${DIRECTORY} ${L} ${H_nu} ${H_e} ${H_eta1}
done

python3 New_Bootstrap_Energy.py ${BASEDIR} ${H_blow} ${H_bhigh} ${nbeta} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho}_alpha${H_alpha} ${nMAX} ${H_init} ${flag} ${LList[@]}
python3 New_Bootstrap_HelicityModulus.py ${BASEDIR} ${H_blow} ${H_bhigh} ${nbeta} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho}_alpha${H_alpha} ${nMAX}  ${H_init} ${flag} ${LList[@]}
python3 New_Bootstrap_Magnetization.py ${BASEDIR} ${H_blow} ${H_bhigh} ${nbeta} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho}_alpha${H_alpha} ${nMAX} ${H_init} ${flag} ${LList[@]}
#python3 New_Bootstrap_DualStiffness.py ${BASEDIR} ${H_blow} ${H_bhigh} ${nbeta} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho}_alpha${H_alpha} ${nMAX} ${H_init} ${flag} ${LList[@]}

#python3 New_Bootstrap_Magnetization_phase.py ${BASEDIR} ${H_blow} ${H_bhigh} ${nbeta} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${nMAX} ${H_init} ${flag} ${LList[@]}
