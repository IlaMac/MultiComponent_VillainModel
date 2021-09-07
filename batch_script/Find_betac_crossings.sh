#!/bin/bash


############# Parameters of the Hamiltonian ##################
H_rho=1
#H_alpha=1
H_eta1=0
H_eta2=0.1
H_e=0
H_h=1
H_nu=0.4
H_blow=0.375
H_bhigh=0.379
bmin=0.376
bmax=0.3785
nMAX=30
H_init="1"


LList=("8 10 12 16 20 24 32 40 44")


BASEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_${H_e}/nu_${H_nu}/h_${H_h}"
#BASEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_${H_e}/nu_${H_nu}/eta2_${H_eta2}/h_${H_h}"


python3 Find_betac_crossings.py ${BASEDIR} ${H_blow} ${H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${nMAX} ${H_init} ${bmin} ${bmax} ${LList[@]}
