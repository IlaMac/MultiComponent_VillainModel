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
            fscanf(fin, "%lf" , &Hp.rho);
            fscanf(fin, "%lf" , &Hp.alpha);
            fscanf(fin, "%lf" , &Hp.rho2);
            fscanf(fin, "%lf" , &Hp.eta1);
            fscanf(fin, "%lf" , &Hp.eta2);
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
        Hp.alpha=1;
        /* To have finite amplitude gradients within the system*/
        Hp.rho2=0;
        Hp.eta1=0;
        Hp.eta2=100;
        Hp.e=0;
        Hp.h= 1.0;
        Hp.nu=0;
        Hp.b_low=0.1;
        Hp.b_high=0.35;
        Hp.init=1;
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
	        //fscanf(fin, "%lf", &MCp.lbox_theta);
            fscanf(fin, "%lf", &MCp.lbox_A);
            fscanf(fin, "%d", &MCp.lbox);
            fscanf(fin, "%d", &MCp.lbox_coupled);
            fscanf(fin, "%d", &MCp.nMAX);
            fclose(fin);
        }
    }else{
        MCp.nmisu=1000;
        MCp.tau=1;
        MCp.n_autosave=2000; //not used now
        //MCp.lbox_theta=C_PI;
        MCp.lbox_A=0.1;
        MCp.lbox=10;
        MCp.lbox_coupled=10;
        MCp.nMAX=30;
    }


}

void initialize_lattice(const std::vector<Node> &Site, const fs::path & directory_read, int RESTART, struct H_parameters &Hp){

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
            for(auto & s : Site){
                fread(s.Psi.data(), sizeof(int), NC, fPsi);
                fread(s.A.data(), sizeof(double), DIM, fA);
            }
            fclose(fA);
            fclose(fPsi);
        }
    }else{
        if(Hp.init==0) {
            for(auto & s : Site){
                for(auto & p : s.Psi){
                    p = 0;
                }
                for(auto & a : s.A){
                    a = 0;
                }
            }
        }
        else if(Hp.init==1) {
            for(auto & s : Site){
                for(auto & p : s.Psi){
                    p = rn::uniform_integer_box(0, MaxP-1) - 0.5*(MaxP-1);
                }
                for(auto & a : s.A){
                    a = 0;
                }
            }
        }
        else if(Hp.init==2) {
            /*This initial conditions correspond to the case where the phase difference is ordered and the phase sum is not (quartic metal) */
            for(auto & s : Site){
                s.Psi[0] = rn::uniform_integer_box(0, MaxP-1) - 0.5*(MaxP-1);
                s.Psi[1]=s.Psi[0] + (MaxP)*0.25;  
                
                for(auto & a : s.A){
                    a = 0;
                }              
            }
        }
        else if(Hp.init==3) {
            /* This initial conditions correspond to the case where the phase of the component 1 is ordered,
             * while both the phase difference and the phase of the component 2 are disordered  */
            for(auto & s : Site){
                double shift_phase=rn::uniform_integer_box(0, MaxP-1) - 0.5*(MaxP-1);
                s.Psi[0] = 0;
                s.Psi[1]=s.Psi[0] + shift_phase;
                for(auto & a : s.A){
                    a = 0;
                }              
            }
        }else if((Hp.init>3) || (Hp.init<0) ) {
            /*This initial conditions correspond to the case where both the phase difference and the phase sum are ordered */
            for(auto & s : Site){
                s.Psi[0] = 0;
                s.Psi[1]=s.Psi[0] + (MaxP)*0.25;

                for(auto & a : s.A){
                    a = 0;
                }
            }
        }
    }
}




