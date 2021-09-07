#!/bin/bash

################# This script computes for a set of L (the array is defined in the .py files): ####################################
#	-the maximum transient time among all the observables
#	-the maximum autocorrelation time among all the observables
#	-Given these two times it perform a bootstrap resampling to compute the mean value and the variance of the observables 
###################################################################################################################################



############# Parameters of the Hamiltonian ##################
H_rho=1
#H_alpha=1
H_eta1=0
H_eta2=0.1
H_e=0
H_h=1
H_nu=0.5
H_blow=0.398
H_bhigh=0.404
nMAX=30
H_init="1"
flag=1 #vecchia nomenclatura cartella senza _init${H_init}

nbeta=64

#LList=("8 10 12 16 20 24 32")
LList=("40 44 64")

#BASEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_${H_e}/nu_${H_nu}/h_${H_h}"
BASEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_${H_e}/nu_${H_nu}/eta2_${H_eta2}/h_${H_h}"


for L in $LList; do
    old_dirname=$BASEDIR/L${L}_rho${H_rho}_eta1${H_eta1}_eta2${H_eta2}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_nMAX${nMAX}_init${H_init}
    new_dirname=$BASEDIR/L${L}_rho${H_rho}_alpha1_eta1${H_eta1}_eta2${H_eta2}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_nMAX${nMAX}_init${H_init}
    mv ${old_dirname} ${new_dirname}
done