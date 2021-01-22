//
// Created by ilaria on 2020-12-17.
//

#include "villain_MC.h"


void metropolis_villain(struct Node* Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta, struct Villain &vil){

    int ix, iy, iz;
    int ip, im, alpha, vec, i;
    int n_var, start=0.5*(MaxP*MaxP-1);
    int new_int_phase[NC]={0};
    //int old_int_phase;
    int arg_F_new[NC][3]={0};
    int arg_B_new[NC][3]={0}; /*Forward and Backward updated phases*/
    int arg_F_old[NC][3]={0};
    int arg_B_old[NC][3]={0}; /*Forward and Backward updated phases*/
    double acc_rate=0.5, acc_theta=0., rand;
    double dE;

    class_tic_toc t_localHtheta(true,5,"local_Htheta");


    for (iz= 0; iz < Lz; iz++) {
        for (iy = 0; iy < Ly; iy++) {
            for (ix = 0; ix < Lx; ix++) {
                i = ix + Lx * (iy + iz * Ly);

                /*******PHASE ONLY UPDATE**************/
                for (alpha = 0; alpha < NC; alpha++) {
                    n_var = rn::uniform_integer_box(1, MCp.lbox);
                    //std::cout<< n_var<< " box: "<< MCp.lbox<< std::endl;
                    rand = rn::uniform_real_box(0, 1);
                    if(rand<0.5){
                        new_int_phase[alpha]= arg(Site[i].Psi[alpha] + n_var, MaxP);}
                    else{
                        new_int_phase[alpha]= arg(Site[i].Psi[alpha] - n_var, MaxP);
                    }
                    //std::cout<< "Component: "<< alpha << " Old phase: "<<  Site[i].Psi[alpha] << " New phase: " << new_int_phase[alpha] << " increment: " << n_var << std::endl;
                        for (vec = 0; vec < 3; vec++) {
                            if (vec == 0) {
                                ip = mod(ix + 1, (int)Lx) + Lx * (iy + iz * Ly);
                                im = mod(ix - 1, (int)Lx) + Lx * (iy + iz * Ly);
                            }
                            if (vec == 1) {
                                ip = ix + Lx * (mod(iy + 1, (int)Ly) + iz * Ly);
                                im = ix + Lx * (mod(iy - 1, (int)Ly) + iz * Ly);
                            }
                            if (vec == 2) {
                                ip = ix + Lx * (iy + mod(iz + 1, (int)Lz) * Ly);
                                im = ix + Lx * (iy + mod(iz - 1, (int)Lz) * Ly);
                            }
                            arg_F_new[alpha][vec] = arg( (Site[ip].Psi[alpha] - new_int_phase[alpha]), MaxP);
                            arg_B_new[alpha][vec] = arg( (new_int_phase[alpha] - Site[im].Psi[alpha]), MaxP);
                            arg_F_old[alpha][vec] = arg( (Site[ip].Psi[alpha] - Site[i].Psi[alpha]), MaxP);
                            arg_B_old[alpha][vec] = arg( (Site[i].Psi[alpha] - Site[im].Psi[alpha]), MaxP);

                           // std::cout<< "arg F new without arg : "<< Site[ip].Psi[alpha] - new_int_phase[alpha]  << " arg F old without arg: "<<  Site[ip].Psi[alpha] - Site[i].Psi[alpha]  << std::endl;

                        }
                    //std::cout<< "BEFORE B_old[0][0]: "<< arg_B_old[0][0]<< " B_old[0][1]: "<< arg_B_old[0][1]<< " B_old[0][2]: "<< arg_B_old[0][2]  << std::endl;
                    //std::cout<< "BEFORE B_new[0][0]: "<< arg_B_new[0][0]<< " B_new[0][1]: "<< arg_B_new[0][1]<< " B_new[0][2]: "<< arg_B_new[0][2]  << std::endl;
                }
                    /*Specific for the two component case*/
                    /*Phase updating phase componet 1*/
                    dE=0;
                    for (vec = 0; vec < 3; vec++) {
                        dE += (vil.potential[start + arg_B_new[0][vec] + MaxP * arg_B_old[1][vec]]
                              +vil.potential[start + arg_F_new[0][vec] + MaxP * arg_F_old[1][vec]]
                              -vil.potential[start + arg_B_old[0][vec] + MaxP * arg_B_old[1][vec]]
                              -vil.potential[start + arg_F_old[0][vec] + MaxP * arg_F_old[1][vec]]);
                    }
                    rand = rn::uniform_real_box(0, 1);
                    //
                    //Boltzmann weight: exp(-\beta dH) dH= 1/beta dE
                    if (rand <exp(-my_beta*dE)) {
                        Site[i].Psi[0] = new_int_phase[0];
                        acc_theta++;

                        for(vec=0; vec<3;vec++){
                            arg_B_old[0][vec]=arg_B_new[0][vec];
                            arg_F_old[0][vec]=arg_F_new[0][vec];
                        }
                    }

                /*Phase updating phase componet 2*/
//               dE=0;
//               for (vec = 0; vec < 3; vec++) {
//                   dE += (vil.potential[start + arg_B_old[0][vec] + MaxP * arg_B_new[1][vec]]
//                         +vil.potential[start + arg_F_old[0][vec] + MaxP * arg_F_new[1][vec]]
//                         -vil.potential[start + arg_B_old[0][vec] + MaxP * arg_B_old[1][vec]]
//                         -vil.potential[start + arg_F_old[0][vec] + MaxP * arg_F_old[1][vec]]);
//               }
//                rand = rn::uniform_real_box(0, 1);
//                //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
//                if (rand < exp(-my_beta*dE)) {
//                    Site[i].Psi[1] = new_int_phase[1];
//                    acc_theta++;
//                }
            }
        }
    }

    acc_theta=(double) acc_theta/N; /*to add NC later*/
    //std::cout << acc_theta/acc_rate << "beta: "<<my_beta<<" lbox: "<< MCp.lbox <<std::endl;
    if(acc_theta/acc_rate>1){MCp.lbox+=1;}
    if(acc_theta/acc_rate<1){MCp.lbox-=1;}
    //if(MCp.lbox>(0.25*(MaxP-1))){MCp.lbox=0.25*(MaxP-1);}

}


/*arg_phase leads x between 0 and 2pi*/
double arg_phase(double x){
    while(x < 0 ){
        x+=2*M_PI;
    }
    while(x>=2*M_PI){
        x-=2*M_PI;
    }
    return x;
}

int arg(int x, int Max){

    while(x < -(Max -1)/2){
        x+=(Max-1);
    }
    while(x > (Max -1)/2){
        x-=(Max-1);
    }
    return x;
}