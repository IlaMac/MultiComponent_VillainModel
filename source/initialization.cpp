//
// Created by ilaria on 2019-11-13.
//
#include "initialization.h"
#include "robust_filesystem.h"

void initialize_Hparameters(struct H_parameters &Hp, const fs::path & directory_parameters){

    fs::path hp_init_file = directory_parameters / "HP_init.txt";
    if(fs::exists(hp_init_file)){
        FILE *fin= nullptr;
        if((fin=fopen(hp_init_file.c_str(), "r"))) {
            fscanf(fin, "%d" , &Hp.rho);
            fscanf(fin, "%d" , &Hp.eta);
            fscanf(fin, "%lf" , &Hp.e);
            fscanf(fin, "%lf" , &Hp.h);
            fscanf(fin, "%lf" , &Hp.nu);
            fscanf(fin, "%lf" , &Hp.b_low);
            fscanf(fin, "%lf" , &Hp.b_high);
            fscanf(fin, "%d" , &Hp.init);
            fclose(fin);
            //With this modification Hp.beta in not anymore part of the Hamiltonian parameters list
        }
    }else{
        Hp.rho=1;
        Hp.eta=0;
        Hp.e=0;
        Hp.h= 1.0;
        Hp.nu=0.1;
        Hp.b_low=0.45;
        Hp.b_high=0.46;
        Hp.init=0;
    }

}

void initialize_MCparameters(struct MC_parameters &MCp, const fs::path & directory_parameters){

    fs::path mc_init_file = directory_parameters / "MC_init.txt";
    if(fs::exists(mc_init_file)){
        FILE *fin= nullptr;
        if((fin=fopen(mc_init_file.c_str(), "r"))) {
            fscanf(fin, "%d", &MCp.nmisu);
            fscanf(fin, "%d", &MCp.tau);
            fscanf(fin, "%d", &MCp.n_autosave);
	        fscanf(fin, "%lf", &MCp.lbox_theta);
            fscanf(fin, "%d", &MCp.lbox);
            fscanf(fin, "%d", &MCp.nMAX);
            fclose(fin);
        }
    }else{
        MCp.nmisu=100;
        MCp.tau=1;
        MCp.n_autosave=2000; //not used now
        MCp.lbox_theta=C_PI;
        //MCp.lbox_A=0.1;
        MCp.lbox=10;
        MCp.nMAX=30;
    }


}

void initialize_lattice(struct Node* Site, const fs::path & directory_read, int RESTART, struct H_parameters &Hp){

    unsigned int i, alpha;
    fs::path psi_init_file = directory_read / "Psi_restart.bin";
    fs::path a_init_file = directory_read / "A_restart.bin";

    if(RESTART==1){
        fs::path psi_init_file = directory_read / "Psi_restart.bin";
        fs::path a_init_file = directory_read / "A_restart.bin";
    }
    else if (RESTART==2){
        fs::path psi_init_file = directory_read / "Psi_final.bin";
        fs::path a_init_file = directory_read / "A_final.bin";
    }

    if((fs::exists(psi_init_file)) and (fs::exists(a_init_file)) and (RESTART!=0)){

        FILE *fPsi= nullptr;
        FILE *fA= nullptr;
        if((fPsi=fopen(psi_init_file.c_str(), "r")) and (fA=fopen(a_init_file.c_str(), "r")) ) {
            for (i = 0; i < N; i++) {
                fread(Site[i].Psi, sizeof(struct O2), NC, fPsi);
                fread(Site[i].A, sizeof(double), 3, fA);
            }
            fclose(fA);
            fclose(fPsi);
        }
    }else{
        if(Hp.init==0) {
            for (i = 0; i < N; i++) {
                for (alpha = 0; alpha < NC; alpha++) {
                    Site[i].Psi[alpha].r =1;
                    Site[i].Psi[alpha].t = 0.;
                    polar_to_cartesian(Site[i].Psi[alpha]);
                }
            }
        }
        else if(Hp.init!=0) {
            for (i = 0; i < N; i++) {
                for (alpha = 0; alpha < NC; alpha++) {
                    Site[i].Psi[alpha].r = 1;
                    Site[i].Psi[alpha].t = rn::uniform_real_box(0, C_TWO_PI);
                    polar_to_cartesian(Site[i].Psi[alpha]);
                }
            }
        }

    }

}




