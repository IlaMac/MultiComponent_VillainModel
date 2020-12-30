//
// Created by ilaria on 2020-12-17.
//

#include "villain_MC.h"


void metropolis_villain(struct Node* Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta, struct Villain &vil){

    unsigned int ix, iy, iz, ip, im, alpha, beta, vec, i;
    double dp=C_TWO_PI/MaxP, inv_dp=MaxP/C_TWO_PI;
    int n_var, start=0.5*(MaxP*MaxP-1);
    int new_int_phase[NC]={0};
    int old_int_phase;
    double new_phase[NC]={0};
    //struct Node* New_Site;
    //struct Node* Old_Site;
    int arg_F_new[NC][3]={0};
    int arg_B_new[NC][3]={0}; /*Forward and Backward updated phases*/
    int arg_F_old[NC][3]={0};
    int arg_B_old[NC][3]={0}; /*Forward and Backward updated phases*/
    double acc_rate=0.5, acc_theta=0., rand;
    double newE, oldE, dE;

    class_tic_toc t_localHtheta(true,5,"local_Htheta");

    //Old_Site=Site;
    //New_Site=Old_Site;

    for (iz= 0; iz < Lz; iz++) {
        for (iy = 0; iy < Ly; iy++) {
            for (ix = 0; ix < Lx; ix++) {
                i = ix + Lx * (iy + iz * Ly);

                /*******PHASE ONLY UPDATE**************/
                for (alpha = 0; alpha < NC; alpha++) {
                    newE=0.;
                    oldE=0.;
                    memset(new_phase, 0, sizeof(new_phase));
                    n_var = rn::uniform_integer_box(-MCp.lbox, MCp.lbox);
                    old_int_phase= Site[i].Psi[alpha].t*inv_dp;
                    new_int_phase[alpha]= int_arg_phase(old_int_phase + n_var);
                    //new_phase[alpha] = Site[i].Psi[alpha].t + (double) n_var * dp;
                    new_phase[alpha]=dp*new_int_phase[alpha];
                    /*Phases defined in the interval [0, 2pi)*/
                    //new_phase[alpha] = arg_phase(new_phase[alpha]);

                    for (beta = 0; beta < NC; beta++) {

                        for (vec = 0; vec < 3; vec++) {
                            if (vec == 0) {
                                ip = mod(ix + 1, Lx) + Lx * (iy + iz * Ly);
                                im = mod(ix - 1, Lx) + Lx * (iy + iz * Ly);
                            }
                            if (vec == 1) {
                                ip = ix + Lx * (mod(iy + 1, Ly) + iz * Ly);
                                im = ix + Lx * (mod(iy - 1, Ly) + iz * Ly);
                            }
                            if (vec == 2) {
                                ip = ix + Lx * (iy + mod(iz + 1, Lz) * Ly);
                                im = ix + Lx * (iy + mod(iz - 1, Lz) * Ly);
                            }
                            /*arg is an integer between -(MaxP -1)/2 and (MaxP -1)/2*/
                            if(beta==alpha){
                                //new_phase[beta]=new_phase[alpha];
                                new_int_phase[beta]= new_int_phase[alpha];
                            }else{
                                //new_phase[beta]=Site[i].Psi[beta].t;
                                new_int_phase[beta]= Site[i].Psi[beta].t*inv_dp;
                            }
                            arg_F_new[beta][vec] = //ARG( (Site[ip].Psi[beta].t*inv_dp - new_int_phase[beta]), MaxP);
                                    arg( (Site[ip].Psi[beta].t*inv_dp - new_int_phase[beta]));
                            //std::cout<< arg_F_new[beta][vec] << " other way: "<< (arg(Site[ip].Psi[beta].t*inv_dp - new_int_phase[beta])) <<
                            //" before ARG: " << (Site[ip].Psi[beta].t - new_phase[beta])*inv_dp << " new: "<< Site[ip].Psi[beta].t*inv_dp - new_int_phase[beta] << std::endl;
                            arg_B_new[beta][vec] = //ARG( (new_int_phase[beta] - Site[im].Psi[beta].t*inv_dp), MaxP);
                                    arg( (new_int_phase[beta] - Site[im].Psi[beta].t*inv_dp));
                            arg_F_old[beta][vec] = //ARG( (Site[ip].Psi[beta].t - Site[i].Psi[beta].t)*inv_dp, MaxP);
                                    arg( (Site[ip].Psi[beta].t - Site[i].Psi[beta].t)*inv_dp);
                            arg_B_old[beta][vec] = //ARG( (Site[i].Psi[beta].t - Site[im].Psi[beta].t)*inv_dp, MaxP);
                                    arg( (Site[i].Psi[beta].t - Site[im].Psi[beta].t)*inv_dp);
                            }
                    }

                    /*Specific for the two component case*/
                    for (vec = 0; vec < 3; vec++) {
                        //std::cout << start + arg_B_new[0][vec] + MaxP*arg_B_new[1][vec] << " " << start + arg_F_new[0][vec] + MaxP*arg_F_new[1][vec] << " MAX:" << start + (MaxP*(MaxP-1)*0.5 +  (MaxP-1)*0.5) <<std::endl;
                        oldE += vil.potential[start + arg_B_old[0][vec] + MaxP * arg_B_old[1][vec]] +
                                vil.potential[start + arg_F_old[0][vec] + MaxP * arg_F_old[1][vec]];
                        newE += vil.potential[start + arg_B_new[0][vec] + MaxP * arg_B_new[1][vec]] +
                                vil.potential[start + arg_F_new[0][vec] + MaxP * arg_F_new[1][vec]];
                        //std::cout<< "MC Villain Potential Old: "<< vil.potential[start + arg_B_old[0][vec] + MaxP * arg_B_old[1][vec]]<< " arg1: "<< arg_B_old[0][vec]  << " arg2: "<< arg_B_old[1][vec]<< std::endl;
                        //std::cout<< "MC Villain Potential New: "<< vil.potential[start + arg_B_new[0][vec] + MaxP * arg_B_new[1][vec]]<< " arg1: "<< arg_B_new[0][vec]  << " arg2: "<< arg_B_new[1][vec]<< std::endl;

                    }
                    dE = newE - oldE;

                    if (dE < 0) {
                        Site[i].Psi[alpha].t = new_phase[alpha];
                        acc_theta++;
//                        //mis.E+=dE;
                    } else {
                        rand = rn::uniform_real_box(0, 1);
                        //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
                        if ((rand) < exp(-my_beta*dE)) {
                            Site[i].Psi[alpha].t = new_phase[alpha];
                            acc_theta++;
                            //mis.E+=dE;
                        }
                    }
                }
            }
        }
    }

    acc_theta=(double) acc_theta/(NC*N);
    //std::cout << acc_theta/acc_rate << "beta: "<<my_beta<<" lbox: "<< MCp.lbox <<std::endl;
    if(acc_theta/acc_rate>1){MCp.lbox+=1;}
    if(acc_theta/acc_rate<1){MCp.lbox-=1;}
    if(MCp.lbox>(0.25*(MaxP-1))){MCp.lbox=0.25*(MaxP-1);}

}


/*arg_phase leads x between 0 and 2pi*/
double arg_phase(double x){
    while(x < 0 ){
        x+=C_TWO_PI;
    }
    while(x>=C_TWO_PI){
        x-=C_TWO_PI;
    }
    return x;
}
int int_arg_phase(int x){
    while(x < 0 ){
        x+=MaxP-1;
    }
    while(x>MaxP-1){
        x-=MaxP-1;
    }
    return x;
}
int arg(int x){

    while(2*x < -(MaxP -1)){
        x+=(MaxP-1);
    }
    while(2*x > (MaxP -1)){
        x-=(MaxP-1);
    }
    return x;
}