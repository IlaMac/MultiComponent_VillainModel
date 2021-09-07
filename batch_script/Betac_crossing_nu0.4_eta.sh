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


eta2_list=("0.1 1 5 10 15")

HOMEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/"

rm -rf ${HOMEDIR}/betac_z2_e${H_e}_nu${H_nu}_vs_eta.txt
rm -rf ${HOMEDIR}/betac_u1_e${H_e}_nu${H_nu}_vs_eta.txt

for H_eta2 in ${eta2_list}; do
    
    var=${H_eta2}
    echo $var

    if [ "${H_eta2}" == "0.1" ]; then 

        Z2_H_blow=0.376
        Z2_H_bhigh=0.379
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.376
        Z2_bmax=0.379
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}
    fi


    if [ "${H_eta2}" == "0.5" ]; then 


        Z2_H_blow=0.36
        Z2_H_bhigh=0.376
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.36
        Z2_bmax=0.376
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}

    fi


    if [ "${H_eta2}" == "1" ]; then 

        Z2_H_blow=0.365
        Z2_H_bhigh=0.367
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.365
        Z2_bmax=0.367
        U1_bmin=0.365
        U1_bmax=${Z2_bmax}

    fi

    if [ "${H_eta2}" == "5" ]; then 

        Z2_H_blow=0.32
        Z2_H_bhigh=0.345
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.32
        Z2_bmax=0.345
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}

    fi

    if [ "${H_eta2}" == "10" ]; then 

	H_init="3"

        Z2_H_blow=0.32
        Z2_H_bhigh=0.33
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.32
        Z2_bmax=0.33
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}

    fi

    if [ "${H_eta2}" == "15" ]; then 

        H_init="3"

        LList=("8 10 12 16 20")

        Z2_H_blow=0.317
        Z2_H_bhigh=0.327
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.317
        Z2_bmax=0.327
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}

    fi

        if [ "${H_eta2}" == "20" ]; then

        H_init="3"

	LList=("8 10 12")

        Z2_H_blow=0.29
        Z2_H_bhigh=0.34
        U1_H_blow=${Z2_H_blow}
        U1_H_bhigh=${Z2_H_bhigh}
        Z2_bmin=0.3
        Z2_bmax=0.33
        U1_bmin=${Z2_bmin}
        U1_bmax=${Z2_bmax}

    fi


    BASEDIR="/home/ilaria/Desktop/MultiComponent_VillainModel/Output_Villain_2C/Model_Sym/e_${H_e}/nu_${H_nu}/eta2_${H_eta2}/h_${H_h}"


    python3 Find_Z2betac_crossings.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]}>> ${HOMEDIR}/betac_z2_e${H_e}_nu${H_nu}_vs_eta.txt
    python3 Find_U1betac_crossings.py ${BASEDIR} ${U1_H_blow} ${U1_H_bhigh} ${H_e} ${H_h} ${H_nu} ${H_eta1} ${H_eta2} ${H_rho} ${H_alpha} ${nMAX} ${H_init} ${U1_bmin} ${U1_bmax} ${var} ${LList[@]}>> ${HOMEDIR}/betac_u1_e${H_e}_nu${H_nu}_vs_eta.txt

done
