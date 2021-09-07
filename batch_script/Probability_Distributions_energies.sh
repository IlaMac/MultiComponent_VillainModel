#!/bin/bash


############# Parameters of the Hamiltonian ##################
H_rho=1
H_alpha=1
H_e=0
H_eta1=0
H_eta2=10
H_h=1
H_nu=0.4
nMAX=30
H_init="1"
H_blow=0.42
H_bhigh=0.44
nbeta=64

LList=("32")

HOMEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/"
BASEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_${H_e}/nu_${H_nu}/eta2_${H_eta2}/h_${H_h}"

python3 Probability_Distributions_energies.py ${BASEDIR} ${H_blow} ${H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${nbeta} ${LList[@]} 
