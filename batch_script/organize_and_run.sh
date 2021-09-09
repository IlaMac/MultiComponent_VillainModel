#!/bin/bash

BASEDIR=${HOME}/MultiComponent_VillainModel
SCRIPT_DIR=${BASEDIR}/batch_script

cd /tmp/

if [ ! -d ./SOutput_x_ilaria ]; then
   mkdir -p Output_x_ilaria
fi

#RESTART=0-> Restart from scratch
#RESTART=1-> Restart from interrupted run
#RESTART=2-> Restart from the previois final scenario

RESTART=0

time_limit="7-00:00:00"

LLIST="8 10 12 16 20 24 32"
#LLIST="8 10 12 16 20 24"
############ Parameters of the Hamiltonian ---> HP_init.txt in a directory whose name contains the main parameters values##################
H_rho=1
H_alpha=1
H_eta1=0
H_eta2=100
H_e=0
H_h=1
H_nu=0
H_blow=0.05
H_bhigh=0.3
H_init=4
#If H_init=0: phases initialized to zero -> phase sum ordered, phase difference not in its minimum value; 
#H_init=1: phases initialized randomly; 
#H_init=2: the phase difference is ordered and the phase sum is not (quartic metal);
#H_init=3: the phase of component1 is ordered while that of component 2 is not;
#H_init=4: both phase difference and phase sum are ordered

############ Parameters for the Monte Carlo simulations --> MC_init.txt#####################

Nmisu=300000
ntau=32
nautosave=100000
theta_box=3.141592653
#A_box=0.1
l_box=10
l_box_coupled=10
nMAX=30

for L in $LLIST; do

############Creation of the output folder and of the two files of initialization####################

cd ${BASEDIR}/Output_Villain_2C

if [ ! -d ./SModel_Sym ]; then
   mkdir -p Model_Sym
fi

cd Model_Sym

if [ ! -d ./Se_${H_e} ]; then
   mkdir -p e_${H_e}
fi

cd e_${H_e}

if [ ! -d ./Snu_${H_nu} ]; then
   mkdir -p nu_${H_nu}
fi

cd nu_${H_nu}

if [ ! -d ./Seta2_${H_eta2} ]; then
   mkdir -p eta2_${H_eta2}
fi

cd eta2_${H_eta2}

if [ ! -d ./Sh_${H_h} ]; then
   mkdir -p h_${H_h}
fi

cd h_${H_h}

if [ ! -d ./SL${L}_rho${H_rho}_alpha${H_alpha}_eta1${H_eta1}_eta2${H_eta2}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_nMAX${nMAX}_init${H_init} ]; then
   mkdir -p L${L}_rho${H_rho}_alpha${H_alpha}_eta1${H_eta1}_eta2${H_eta2}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_nMAX${nMAX}_init${H_init}
fi

OUTPUT=${BASEDIR}/Output_Villain_2C/Model_Sym/e_${H_e}/nu_${H_nu}/eta2_${H_eta2}/h_${H_h}/L${L}_rho${H_rho}_alpha${H_alpha}_eta1${H_eta1}_eta2${H_eta2}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_nMAX${nMAX}_init${H_init}

cd /tmp/Output_x_ilaria

if [ ! -d ./SModel_Sym ]; then
   mkdir -p Model_Sym
fi

cd Model_Sym

if [ ! -d ./Se_${H_e} ]; then
   mkdir -p e_${H_e}
fi

cd e_${H_e}

if [ ! -d ./Snu_${H_nu} ]; then
   mkdir -p nu_${H_nu}
fi

cd nu_${H_nu}

if [ ! -d ./Seta2_${H_eta2} ]; then
   mkdir -p eta2_${H_eta2}
fi

cd eta2_${H_eta2}

if [ ! -d ./Sh_${H_h} ]; then
   mkdir -p h_${H_h}
fi

cd h_${H_h}

if [ ! -d ./SL${L}_rho${H_rho}_alpha${H_alpha}_eta1${H_eta1}_eta2${H_eta2}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_nMAX${nMAX}_init${H_init} ]; then
   mkdir -p L${L}_rho${H_rho}_alpha${H_alpha}_eta1${H_eta1}_eta2${H_eta2}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_nMAX${nMAX}_init${H_init}
fi

OUTPUT_TEMP=/tmp/Output_x_ilaria/Model_Sym/e_${H_e}/nu_${H_nu}/eta2_${H_eta2}/h_${H_h}/L${L}_rho${H_rho}_alpha${H_alpha}_eta1${H_eta1}_eta2${H_eta2}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_nMAX${nMAX}_init${H_init}

cd ${OUTPUT}

#THE ORDER OF WRITING DOES MATTER
echo $H_rho > HP_init.txt
echo $H_alpha >> HP_init.txt
echo $H_eta1 >> HP_init.txt
echo $H_eta2 >> HP_init.txt
echo $H_e >> HP_init.txt
echo $H_h >> HP_init.txt
echo $H_nu >> HP_init.txt
echo $H_blow >> HP_init.txt
echo $H_bhigh >> HP_init.txt
	echo $H_init >> HP_init.txt

#THE ORDER OF WRITING DOES MATTER
echo $Nmisu > MC_init.txt
echo $ntau >> MC_init.txt
echo $nautosave >> MC_init.txt
echo $theta_box >> MC_init.txt
echo $l_box >> MC_init.txt
echo $l_box_coupled >> MC_init.txt
echo $nMAX >> MC_init.txt

#################Creation of the submit_runs script#########################

jobname="Sym_L${L}_rho${H_rho}_alpha${H_alpha}_eta1${H_eta1}_eta2${H_eta2}_e${H_e}_h${H_h}_nu${H_nu}_bmin${H_blow}_bmax${H_bhigh}_nMAX${nMAX}_init${H_init}"
nnodes=2
ntasks=64 #parallel tempering over ntasks temperatures

#I create ntasks folder: one for each rank.

cd ${OUTPUT}

for ((rank=0; rank<${ntasks}; rank++)); do

if [ ! -d ./Sbeta_${rank} ]; then
   mkdir -p beta_${rank}
fi

done

cd ${OUTPUT_TEMP}

for ((rank=0; rank<${ntasks}; rank++)); do

if [ ! -d ./Sbeta_${rank} ]; then
   mkdir -p beta_${rank}
fi

done

cd ${SCRIPT_DIR}
DIR_PAR="${OUTPUT}"
DIR_PAR_TEMP="${OUTPUT_TEMP}"

#SEED= If I want to repeat exactly a simulation I could initialize the random number generator exactly at the same way

EXECUTE_DIR="../build/Release"

#SBATCH --nodes=${nnodes}               # Number of nodes

echo "#!/bin/bash
#SBATCH --job-name=${jobname}          # Name of the job
#SBATCH --time=${time_limit}               # Allocation time
#SBATCH --mem-per-cpu=2000              # Memory per allocated cpu
#SBATCH --nodes=${nnodes}               # Number of nodes
#SBATCH --ntasks=${ntasks}
#SBATCH --output=${DIR_PAR}/logs/log_${jobname}.o
#SBATCH --error=${DIR_PAR}/logs/log_${jobname}.e

srun ${EXECUTE_DIR}/Villain_2component ${L} ${DIR_PAR} ${DIR_PAR} ${RESTART} &> ${DIR_PAR}/logs/log_${jobname}.o

" >  submit_run

#Submission of the work --> sbatch submit_runs

mkdir -p ${DIR_PAR}/logs

sbatch submit_run

done
