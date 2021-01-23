//
// Created by ilaria on 2020-12-16.
//

#include "villain_potential.h"

void init_villain_potentials(double my_beta, struct Villain &vil,  struct H_parameters &Hp, struct MC_parameters &MCp, const fs::path & directory_write) {

    int n1, n2, arg1, arg2;
    double u1, u2, sum_1, norm, boltz, boltz_H;
    double j1, j2;
    double d1, d2, d11, d12, d22;
    double dp=2*M_PI/MaxP;
    fs::path vpotential_file = directory_write / std::string("Villain_potential.txt");

    FILE *fVPotential= nullptr;
    fVPotential=fopen(vpotential_file.c_str(), "w");

    vil.potential.resize(MaxP*MaxP);
    vil.upotential.resize(MaxP*MaxP);
    vil.d1_potential.resize(MaxP*MaxP);
    vil.d2_potential.resize(MaxP*MaxP);
    vil.d11_potential.resize(MaxP*MaxP);
    vil.d22_potential.resize(MaxP*MaxP);
    vil.d12_potential.resize(MaxP*MaxP);

    for (arg2 = -(MaxP - 1) / 2; arg2 <= (MaxP - 1) / 2; arg2++) {
        for (arg1 = -(MaxP - 1) / 2; arg1 <= (MaxP - 1) / 2; arg1++) {
            norm = 0;
            sum_1 =  0;
            d1=0;
            d2=0;
            d11=0;
            d22=0;
            d12=0;
            for (n2 = -MCp.nMAX; n2 < (MCp.nMAX+1); n2++) {
                for (n1 = -MCp.nMAX; n1 < (MCp.nMAX+1); n1++) {
                    //std::cout<< "u1: "<< u1<< " u2: "<< u2<< " n1: "<< n1 << " n2: "<< n2<< " norm: "<< norm << std::endl;
                    u1=dp*arg1 - 2*M_PI*n1;
                    u2=dp*arg2 - 2*M_PI*n2;
                    //boltz_H=0.5*(Hp.rho*(SQR(dp*arg1-2*M_PI*n1) + SQR(dp*arg2-2*M_PI*n2)) - Hp.nu*SQR(dp*(arg1-arg2)-2*M_PI*(n1-n2)));
                    boltz_H = 0.5*(Hp.rho * ((u1*u1) + (u2*u2)) - Hp.nu*(u1-u2)*(u1-u2));
                    boltz = exp(-my_beta*boltz_H);
                    norm += boltz;
                    sum_1 += boltz_H*boltz;

                    j1=Hp.rho*u1 - Hp.nu*(u1-u2);
                    j2=Hp.rho*u2 - Hp.nu*(u1-u2);

                    d1+= j1*boltz;
                    d2+= j2*boltz;
                    d11+= (Hp.rho - Hp.nu - my_beta*j1*j1)*boltz;
                    d22+= (Hp.rho - Hp.nu - my_beta*j2*j2)*boltz;
                    d12+= (Hp.nu - my_beta*j1*j2)*boltz;

                }
            }
            /*Table Villain potential*/
            vil.potential[OFFSET_POT + arg1 + MaxP*arg2]= -log(norm)/my_beta;
            fprintf(fVPotential, "%d %d %19.12e \n", arg1, arg2, my_beta*vil.potential[OFFSET_POT + arg1 + MaxP*arg2]);
            /*Table for the helicity modulus*/
            vil.d1_potential[OFFSET_POT + arg1 + MaxP*arg2]+= d1/norm;
            vil.d2_potential[OFFSET_POT + arg1 + MaxP*arg2]+= d2/norm;
            vil.d11_potential[OFFSET_POT + arg1 + MaxP*arg2]+=d11/norm + my_beta*(d1/norm)*(d1/norm);
            vil.d22_potential[OFFSET_POT + arg1 + MaxP*arg2]+=d22/norm + my_beta*(d2/norm)*(d2/norm);
            vil.d12_potential[OFFSET_POT + arg1 + MaxP*arg2]+= d12/norm + my_beta*(d1*d2)/(norm*norm);
            /*Calculating beta-derivative of the Villain potential needed for calculating the internal energy (which differs from the energy mean due to the temperature dependence of the Villain model).*/
            vil.upotential[OFFSET_POT + (arg1 + MaxP*arg2)] = sum_1/norm;

        }
    }
    fclose(fVPotential);
}

void init_villainpotential_nnbeta(double beta_np, double beta_nm, struct Villain &vil,  struct H_parameters &Hp, struct MC_parameters &MCp) {

    int n1, n2, arg1, arg2;
    double sum_np, sum_nm, u1, u2;
    double dp=2*M_PI/MaxP;

    vil.potential_bplus.resize(MaxP*MaxP);
    vil.potential_bminus.resize(MaxP*MaxP);
    for (arg2 = -(MaxP - 1) / 2; arg2 <= (MaxP - 1) / 2; arg2++) {
        for (arg1 = -(MaxP - 1) / 2; arg1 <= (MaxP - 1) / 2; arg1++) {
            sum_np = 0;
            sum_nm = 0;
            for (n2 = -MCp.nMAX; n2 < (MCp.nMAX+1); n2++) {
                for (n1 = -MCp.nMAX; n1 < (MCp.nMAX+1); n1++) {
                    u2=dp*arg2 - 2*M_PI*n2;
                    u1=dp*arg1 - 2*M_PI*n1;
                    sum_np+=exp(-0.5*beta_np*(Hp.rho*((u1*u1) + (u2*u2)) - Hp.nu*(u1-u2)*(u1-u2)));
                    sum_nm+=exp(-0.5*beta_nm*(Hp.rho*((u1*u1) + (u2*u2)) - Hp.nu*(u1-u2)*(u1-u2)));
                }
            }
            vil.potential_bplus[OFFSET_POT + arg1 + MaxP*arg2]= -log(sum_np)/beta_np;
            vil.potential_bminus[OFFSET_POT + arg1 + MaxP*arg2]= -log(sum_nm)/beta_nm;
        }
    }
}



