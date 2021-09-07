#!/bin/bash


############# Parameters of the Hamiltonian ##################
H_rho=1
H_alpha=1
H_eta1=0
H_eta2=0.1
H_h=1
H_nu=0
nMAX=30
H_init="1"


LList=("8 10 12 16 20 24 32")


e_list=("0 1 1.5 2 2.5")

HOMEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/"

rm -rf ${HOMEDIR}/betac_z2_eta2${H_eta2}_nu${H_nu}_vs_e.txt
rm -rf ${HOMEDIR}/betac_u1_eta2${H_eta2}_nu${H_nu}_vs_e.txt

for H_e in ${e_list}; do
    
    var=${H_e}
    echo $var
    BASEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_${H_e}/nu_${H_nu}/eta2_${H_eta2}/h_${H_h}"

    if [ "${H_e}" == "0" ]; then 

        Z2_H_blow=0.33
        Z2_H_bhigh=0.34
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.33
        Z2_bmax=0.34
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}

        python3 Find_Z2betac_crossings.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]} >> ${HOMEDIR}/betac_z2_eta2${H_eta2}_nu${H_nu}_vs_e.txt
        python3 Find_U1betac_crossings.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]}>> ${HOMEDIR}/betac_u1_eta2${H_eta2}_nu${H_nu}_vs_e.txt

    fi

       
    if [ "${H_e}" == "0.5" ]; then 
        
        LList=("20 24 32")

        Z2_H_blow=0.35
        Z2_H_bhigh=0.45
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.36
        Z2_bmax=0.38
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}

        python3 Find_Z2betac_crossings.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]} >> ${HOMEDIR}/betac_z2_eta2${H_eta2}_nu${H_nu}_vs_e.txt
        python3 Find_U1betac_crossings_charge.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]}>> ${HOMEDIR}/betac_u1_eta2${H_eta2}_nu${H_nu}_vs_e.txt

    fi


    if [ "${H_e}" == "1" ]; then 

        Z2_H_blow=0.41
        Z2_H_bhigh=0.43
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.41
        Z2_bmax=0.43
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}

        python3 Find_Z2betac_crossings.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]}>> ${HOMEDIR}/betac_z2_eta2${H_eta2}_nu${H_nu}_vs_e.txt
        python3 Find_U1betac_crossings_charge.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]}>> ${HOMEDIR}/betac_u1_eta2${H_eta2}_nu${H_nu}_vs_e.txt

    fi


    if [ "${H_e}" == "1.5" ]; then 

        Z2_H_blow=0.48
        Z2_H_bhigh=0.52
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.48
        Z2_bmax=0.52
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}

        python3 Find_Z2betac_crossings.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]}>> ${HOMEDIR}/betac_z2_eta2${H_eta2}_nu${H_nu}_vs_e.txt
        python3 Find_U1betac_crossings_charge.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]} >> ${HOMEDIR}/betac_u1_eta2${H_eta2}_nu${H_nu}_vs_e.txt

    fi

    if [ "${H_e}" == "2" ]; then 

        Z2_H_blow=0.56
        Z2_H_bhigh=0.6
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.56
        Z2_bmax=0.6
        U1_bmin=0.56
        U1_bmax=0.6

        python3 Find_Z2betac_crossings.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]}>> ${HOMEDIR}/betac_z2_eta2${H_eta2}_nu${H_nu}_vs_e.txt
        python3 Find_U1betac_crossings_charge.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]} >> ${HOMEDIR}/betac_u1_eta2${H_eta2}_nu${H_nu}_vs_e.txt

    fi

    if [ "${H_e}" == "2.5" ]; then 

        Z2_H_blow=0.55
        Z2_H_bhigh=0.65
        U1_H_blow=0.65
        U1_H_bhigh=0.75
        Z2_bmin=0.6
        Z2_bmax=0.65
        U1_bmin=0.7
        U1_bmax=0.75

        python3 Find_Z2betac_crossings.py ${BASEDIR} ${Z2_H_blow} ${Z2_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${Z2_bmin} ${Z2_bmax} ${var} ${LList[@]}>> ${HOMEDIR}/betac_z2_eta2${H_eta2}_nu${H_nu}_vs_e.txt
        python3 Find_U1betac_crossings_charge.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]} >> ${HOMEDIR}/betac_u1_eta2${H_eta2}_nu${H_nu}_vs_e.txt

    fi



done