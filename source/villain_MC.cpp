//
// Created by ilaria on 2020-12-17.
//

#include "villain_MC.h"


void metropolis_villain(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta, struct Villain &vil){

    int ix, iy, iz;
    int ipx, ipy, ipz;
    int ip, im, alpha, vec, i;
    int n_var;
    int new_int_phase[NC]={0};
    int arg_F_new[NC][3]={{0}};
    int arg_B_new[NC][3]={{0}}; /*Forward and Backward updated phases*/
    int arg_F_old[NC][3]={{0}};
    int arg_B_old[NC][3]={{0}}; /*Forward and Backward updated phases*/
    double acc_rate=0.5, acc_theta=0., rand;
    double dE;
    double dE_A;
    double OldA, NewA, dA, acc_A=0.;
    double F_A_new, B_A_new, F_A_old, B_A_old;
    int l, nn_ipl, nn_iml, nn_ipvec_iml;
    double curl2_A_new, curl2_A_old;

    class_tic_toc t_localHtheta(true,5,"local_Htheta");


    for (iz= 0; iz < Lz; iz++) {
        for (iy = 0; iy < Ly; iy++) {
            for (ix = 0; ix < Lx; ix++) {
                i = ix + Lx * (iy + iz * Ly);

                /*******PHASE ONLY UPDATE**************/
                for (alpha = 0; alpha < NC; alpha++) {
                    n_var = rn::uniform_integer_box(1, MCp.lbox);
                    //std::cout<< n_var<< " box: "<< MCp.lbox<< std::endl;
                    //Try to extract a new value of the phase directly in the interval [-(MaxP-1)/2; (MaxP-1/2]
                    rand = rn::uniform_real_box(0, 1);
                    if(rand<0.5){
                        new_int_phase[alpha]= arg(Site[i].Psi[alpha] + n_var, MaxP);}
                    else{
                        new_int_phase[alpha]= arg(Site[i].Psi[alpha] - n_var, MaxP);
                    }
                    //std::cout<< "Component: "<< alpha << " Old phase: "<<  Site[i].Psi[alpha] << " New phase: " << new_int_phase[alpha] << " increment: " << n_var << std::endl;
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
                        //std::cout<< "ix: "<< ix << " iy: "<< iy <<" iz: "<< iz << " vec= "<< vec << " ip: " << ip << " im: "<< im << std::endl;
                         arg_F_new[alpha][vec] = arg( (Site[ip].Psi[alpha] - new_int_phase[alpha]) - inv_dp*Hp.e*Site[i].A[vec], MaxP);
                         arg_B_new[alpha][vec] = arg( (new_int_phase[alpha] - Site[im].Psi[alpha]) - inv_dp*Hp.e*Site[im].A[vec], MaxP);
                         arg_F_old[alpha][vec] = arg( (Site[ip].Psi[alpha] - Site[i].Psi[alpha]) - inv_dp*Hp.e*Site[i].A[vec], MaxP);
                         arg_B_old[alpha][vec] = arg( (Site[i].Psi[alpha] - Site[im].Psi[alpha]) - inv_dp*Hp.e*Site[im].A[vec], MaxP);
                    }
                }
                /******************Specific for the two component case******************/

                /*******************Phase updating phase componet 1*******************/
                dE=0;
                for (vec = 0; vec < DIM; vec++) {
                   dE += (vil.potential.at(OFFSET_POT + arg_B_new[0][vec] + MaxP * arg_B_old[1][vec])
                         +vil.potential.at(OFFSET_POT + arg_F_new[0][vec] + MaxP * arg_F_old[1][vec])
                         -vil.potential.at(OFFSET_POT + arg_B_old[0][vec] + MaxP * arg_B_old[1][vec])
                         -vil.potential.at(OFFSET_POT + arg_F_old[0][vec] + MaxP * arg_F_old[1][vec]));
                }
                rand = rn::uniform_real_box(0, 1);
                dE+= Hp.eta1*cos(dp*(new_int_phase[0] - Site[i].Psi[1])) + Hp.eta2*cos(2*dp*(new_int_phase[0] - Site[i].Psi[1]));
                //Boltzmann weight: exp(-\beta dH) dH= 1/beta dE
                if (rand <exp(-my_beta*dE)) {
                    Site[i].Psi[0] = new_int_phase[0];
                    acc_theta++;
                    for(vec=0; vec<3;vec++){
                        arg_B_old[0][vec]=arg_B_new[0][vec];
                        arg_F_old[0][vec]=arg_F_new[0][vec];
                    }
                }

                /*******************Phase updating phase componet 2*******************/
                dE=0;
                for (vec = 0; vec < DIM; vec++) {
                   dE += (vil.potential.at(OFFSET_POT + arg_B_old[0][vec] + MaxP * arg_B_new[1][vec])
                         +vil.potential.at(OFFSET_POT + arg_F_old[0][vec] + MaxP * arg_F_new[1][vec])
                         -vil.potential.at(OFFSET_POT + arg_B_old[0][vec] + MaxP * arg_B_old[1][vec])
                         -vil.potential.at(OFFSET_POT + arg_F_old[0][vec] + MaxP * arg_F_old[1][vec]));
                }
                dE+= Hp.eta1*cos(dp*(new_int_phase[1] - Site[i].Psi[0])) + Hp.eta2*cos(2*dp*(new_int_phase[1] - Site[i].Psi[0]));
                rand = rn::uniform_real_box(0, 1);
                //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
                if (rand < exp(-my_beta*dE)) {
                    Site[i].Psi[1] = new_int_phase[1];
                    acc_theta++;
                }

            /*******************VECTOR POTENTIAL UPDATE*******************/
                if(Hp.e!=0){
                    for (vec=0; vec<DIM; vec++){
                        OldA= Site[i].A[vec];
                        dA = rn::uniform_real_box(-MCp.lbox_A, MCp.lbox_A);
                        NewA = OldA + dA;
                        ipx=ix;
                        ipy=iy;
                        ipz=iz;
                        if (vec == 0) {
                            ipx=mod(ix + 1, Lx);
                        }
                        if (vec == 1) {
                            ipy=mod(iy + 1, Ly);
                        }
                        if (vec == 2) {
                            ipz = mod(iz + 1, Lz);
                        }
                        ip = ipx + Lx * (ipy + ipz * Ly);
                        for(alpha=0; alpha<NC; alpha++){
                            arg_F_new[alpha][vec] = arg( (Site[ip].Psi[alpha] - Site[i].Psi[alpha]) - inv_dp*Hp.e*NewA, MaxP);
                            arg_F_old[alpha][vec] = arg( (Site[ip].Psi[alpha] - Site[i].Psi[alpha]) - inv_dp*Hp.e*OldA, MaxP);
                        }
                        dE_A= vil.potential.at(OFFSET_POT + arg_F_new[0][vec] + MaxP * arg_F_new[1][vec]) - vil.potential.at(OFFSET_POT + arg_F_old[0][vec] + MaxP * arg_F_old[1][vec]);
                        //curl A
                        curl2_A_new=0.;
                        curl2_A_old=0.;
                        //All the plaquettes involving A_vec(i)
                        for(l=0; l<DIM; l++){
                            if(l!= vec){
                                if(l==0){
                                    nn_ipl= mod(ix + 1, Lx) + Lx * (iy + iz * Ly);
                                    nn_iml= mod(ix - 1, Lx) + Lx * (iy + iz * Ly);
                                    nn_ipvec_iml= mod(ipx - 1, Lx) + Lx * (ipy + ipz * Ly);
                                }else if(l==1){
                                    nn_ipl= ix + Lx * (mod(iy + 1, Ly) + iz * Ly);
                                    nn_iml= ix + Lx * (mod(iy - 1, Ly) + iz * Ly);
                                    nn_ipvec_iml= ipx + Lx * (mod(ipy - 1, Ly) + ipz * Ly);
                                }else if(l==2) {
                                    nn_ipl= ix + Lx * (iy + mod(iz + 1, Lz) * Ly);
                                    nn_iml= ix + Lx * (iy + mod(iz - 1, Lz) * Ly);
                                    nn_ipvec_iml= ipx + Lx * (ipy + mod(ipz - 1, Lz) * Ly);
                                }
                                //plaquette to the right (F="forward") side of newA
                                F_A_new=(NewA + Site[ip].A[l] - Site[nn_ipl].A[vec] - Site[i].A[l]);
                                curl2_A_new+=0.5*(F_A_new*F_A_new);
                                //plaquette to the left (B="backrward") side of newA
                                B_A_new=(Site[nn_iml].A[vec]+ Site[nn_ipvec_iml].A[l] - NewA - Site[nn_iml].A[l]);
                                curl2_A_new+=0.5*(B_A_new*B_A_new);

                                F_A_old=(OldA + Site[ip].A[l] - Site[nn_ipl].A[vec] - Site[i].A[l]);
                                curl2_A_old+=0.5*(F_A_old*F_A_old);
                                B_A_old=(Site[nn_iml].A[vec]+ Site[nn_ipvec_iml].A[l] - OldA - Site[nn_iml].A[l]);
                                curl2_A_old+=0.5*(B_A_old*B_A_old);
                            }
                        }
                        dE_A+= curl2_A_new -curl2_A_old;
                        rand = rn::uniform_real_box(0, 1);
                        //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
                        if (rand < exp(-my_beta*dE_A)) {
                           Site[i].A[vec] = NewA;
                           acc_A++;
                        }
                    }
                }
            }
        }
    }

    acc_theta=(double) acc_theta/(N*NC);
    acc_A=(double) acc_A/(DIM*N);
    //std::cout << acc_theta/acc_rate << "beta: "<<my_beta<<" lbox: "<< MCp.lbox <<std::endl;
    if(acc_theta/acc_rate>1){MCp.lbox+=1;}
    if(acc_theta/acc_rate<1){MCp.lbox-=1;}
    MCp.lbox_A= MCp.lbox_A*((0.5*acc_A/acc_rate)+0.5);

    //if(MCp.lbox>(0.25*(MaxP-1))){MCp.lbox=0.25*(MaxP-1);}

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