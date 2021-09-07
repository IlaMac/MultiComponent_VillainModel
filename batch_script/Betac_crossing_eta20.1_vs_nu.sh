#!/bin/bash


nu_list=("0 0.2 0.4 0.5 0.55 0.6 0.65")
# nu_list=("0.4")
############# Parameters of the Hamiltonian ##################
H_rho=1
H_alpha=1
H_eta1=0
H_eta2=0.1
H_e=0
H_h=1
nMAX=30
H_init="1"


#
LList=("8 10 12 16 20 24 32")
#LList=("8 10 12 16 20 24 32 40 44 64")

HOMEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/"

rm -rf ${HOMEDIR}/betac_z2_vs_nu_eta2${H_eta2}_e${H_e}.txt
rm -rf ${HOMEDIR}/betac_u1_vs_nu_eta2${H_eta2}_e${H_e}.txt

for H_nu in ${nu_list}; do

    var=${H_nu}
    echo ${H_nu}
    LList=("8 10 12 16 20 24 32")

    if [ ${H_nu} == 0 ]; then 
        Z2_H_blow=0.33
        Z2_H_bhigh=0.34
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.33
        Z2_bmax=0.34
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}

    fi

    if [ "${H_nu}" == "0.2" ]; then 
        Z2_H_blow=0.34
        Z2_H_bhigh=0.35
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.34
        Z2_bmax=0.35
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}
    fi

    if [ "${H_nu}" == "0.4" ]; then 
        Z2_H_blow=0.376
        Z2_H_bhigh=0.379
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}       
        Z2_bmin=0.376
        Z2_bmax=0.379
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}
    fi

      
    if [ "${H_nu}" == "0.5" ]; then 
        Z2_H_blow=0.398
        Z2_H_bhigh=0.404
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.398
        Z2_bmax=0.404
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}
    fi

    if [ "${H_nu}" == "0.55" ]; then 
        Z2_H_blow=0.406
        Z2_H_bhigh=0.42
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.406
        Z2_bmax=0.42
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}
    fi

    if [ "${H_nu}" == "0.6" ]; then 
        Z2_H_blow=0.402
        Z2_H_bhigh=0.414
        Z2_bmin=0.402
        Z2_bmax=0.414

        U1_H_blow=0.43
        U1_H_bhigh=0.45
        U1_bmin=0.43
        U1_bmax=0.45
    fi

    if [ "${H_nu}" == "0.65" ]; then 
        Z2_H_blow=0.392
        Z2_H_bhigh=0.402
        Z2_bmin=0.392
        Z2_bmax=0.402

        U1_H_blow=0.475
        U1_H_bhigh=0.5
        U1_bmin=0.475
        U1_bmax=0.5
    fi

    BASEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_${H_e}/nu_${H_nu}/eta2_${H_eta2}/h_${H_h}"

    #echo ${BASEDIR} ${Z2_H_blow} ${Z2_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${nMAX} ${H_init} ${Z2_bmin} ${Z2_bmax} ${var} ${LList[@]}
    python3 Find_Z2betac_crossings.py ${BASEDIR} ${Z2_H_blow} ${Z2_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${Z2_bmin} ${Z2_bmax} ${var} ${LList[@]} >> ${HOMEDIR}/betac_z2_vs_nu_eta2${H_eta2}_e${H_e}.txt
    python3 Find_U1betac_crossings.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]} >> ${HOMEDIR}/betac_u1_vs_nu_eta2${H_eta2}_e${H_e}.txt

done
