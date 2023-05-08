#!/bin/bash


############# Parameters of the Hamiltonian ##################
H_rho=1
H_alpha=1
H_eta1=0
H_e=0
H_h=1
H_nu=0.4
nMAX=30
H_init="1"


LList=("8 10 12 16 20 24 32")


eta2_list=("0.1 10 20 50 100")

HOMEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/"

rm -rf ${HOMEDIR}/betac_z2_e${H_e}_nu${H_nu}_vs_eta.txt
rm -rf ${HOMEDIR}/betac_u1_e${H_e}_nu${H_nu}_vs_eta.txt

for H_eta2 in ${eta2_list}; do
    
    var=${H_eta2}
    echo $var

    if [ "${H_eta2}" == "0.1" ]; then 

        Z2_H_blow=0.3755
        Z2_H_bhigh=0.378
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.3755
        Z2_bmax=0.378
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}
    fi


    if [ "${H_eta2}" == "10" ]; then 

        LList=("10 12 16 20 24 32")

        Z2_H_blow=0.32
        Z2_H_bhigh=0.34
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.32
        Z2_bmax=0.34
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}

    fi


    if [ "${H_eta2}" == "20" ]; then 

        LList=("10 12 16 20 24 32")

        Z2_H_blow=0.315
        Z2_H_bhigh=0.335
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.315
        Z2_bmax=0.335
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}

    fi

    if [ "${H_eta2}" == "50" ]; then 

        LList=("10 12 16 20 24 32")

        Z2_H_blow=0.3
        Z2_H_bhigh=0.33
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.31
        Z2_bmax=0.32
        U1_bmin=0.315
        U1_bmax=0.33

    fi

    if [ "${H_eta2}" == "100" ]; then 

        LList=("10 12 16 20 24 32")

        Z2_H_blow=0.3
        Z2_H_bhigh=0.33
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.31
        Z2_bmax=0.32
        U1_bmin=0.315
        U1_bmax=0.33

    fi



    BASEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_${H_e}/nu_${H_nu}/eta2_${H_eta2}/h_${H_h}"


    python3 Find_Z2betac_crossings.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]} >> ${HOMEDIR}/betac_z2_e${H_e}_nu${H_nu}_vs_eta.txt
    python3 Find_U1betac_crossings.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]} >> ${HOMEDIR}/betac_u1_e${H_e}_nu${H_nu}_vs_eta.txt

done
