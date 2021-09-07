#!/bin/bash


alpha_list=("0.7 0.8 0.9 1")

############# Parameters of the Hamiltonian ##################
H_rho=1
H_eta1=0
H_eta2=0.1
H_e=0
H_nu=0.6
H_h=1
nMAX=30
H_init="1"


LList=("8 10 12 16 20 24 32")

HOMEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/"

rm -rf ${HOMEDIR}/betac_z2_vs_alpha_nu${H_nu}_eta2${H_eta2}_e${H_e}.txt
rm -rf ${HOMEDIR}/betac_u1_vs_alpha_nu${H_nu}_eta2${H_eta2}_e${H_e}.txt

for H_alpha in ${alpha_list}; do

    var=${H_alpha}
    echo ${H_alpha}

    if [ "${H_alpha}" == "1" ]; then 
        Z2_H_blow=0.4
        Z2_H_bhigh=0.41
        Z2_bmin=0.4
        Z2_bmax=0.41

        U1_H_blow=0.43
        U1_H_bhigh=0.45
        U1_bmin=0.43
        U1_bmax=0.45
    fi

    if [ "${H_alpha}" == "0.9" ]; then 

        Z2_H_blow=0.4
        Z2_H_bhigh=0.44
        Z2_bmin=0.4
        Z2_bmax=0.44

        U1_H_blow=0.45
        U1_H_bhigh=0.525
        U1_bmin=0.45
        U1_bmax=0.525
    fi

    if [ "${H_alpha}" == "0.8" ]; then 

        Z2_H_blow=0.4
        Z2_H_bhigh=0.5
        Z2_bmin=0.4
        Z2_bmax=0.5

        U1_H_blow=0.5
        U1_H_bhigh=0.65
        U1_bmin=0.5
        U1_bmax=0.65
    fi

    if [ "${H_alpha}" == "0.7" ]; then 

        Z2_H_blow=0.4
        Z2_H_bhigh=0.5
        Z2_bmin=0.4
        Z2_bmax=0.5

        U1_H_blow=0.625
        U1_H_bhigh=0.7
        U1_bmin=0.625
        U1_bmax=0.7
    fi


    BASEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_${H_e}/nu_${H_nu}/eta2_${H_eta2}/h_${H_h}"

    #echo ${BASEDIR} ${Z2_H_blow} ${Z2_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${nMAX} ${H_init} ${Z2_bmin} ${Z2_bmax} ${var} ${LList[@]}
    python3 Find_Z2betac_crossings.py ${BASEDIR} ${Z2_H_blow} ${Z2_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${Z2_bmin} ${Z2_bmax} ${var} ${LList[@]} >> ${HOMEDIR}/betac_z2_vs_alpha_nu${H_nu}_eta2${H_eta2}_e${H_e}.txt
    python3 Find_U1betac_crossings.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]} >> ${HOMEDIR}/betac_u1_vs_alpha_nu${H_nu}_eta2${H_eta2}_e${H_e}.txt

done
