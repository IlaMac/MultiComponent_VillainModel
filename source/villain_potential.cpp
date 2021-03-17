//
// Created by ilaria on 2020-12-16.
//

#include "villain_potential.h"

void init_villain_potentials(double my_beta, double beta_np, double beta_nm, struct Villain &vil,  struct H_parameters &Hp, struct MC_parameters &MCp, const fs::path & directory_write) {

    double u1, u2, sum_1, norm, boltz, boltz_H;
    double j1, j2;

    fs::path vpotential_file = directory_write / std::string("Villain_potential.txt");

    FILE *fVPotential;
    fVPotential=fopen(vpotential_file.c_str(), "w");

    vil.potential.resize(MaxP*MaxP);
    vil.upotential.resize(MaxP*MaxP);
    vil.d1_potential.resize(MaxP*MaxP);
    vil.d2_potential.resize(MaxP*MaxP);
    vil.d11_potential.resize(MaxP*MaxP);
    vil.d22_potential.resize(MaxP*MaxP);
    vil.d12_potential.resize(MaxP*MaxP);
    vil.potential_bplus.resize(MaxP*MaxP);
    vil.potential_bminus.resize(MaxP*MaxP);

    for (int arg2 = -(MaxP - 1) / 2; arg2 <= (MaxP - 1) / 2; arg2++) {
        for (int arg1 = -(MaxP - 1) / 2; arg1 <= (MaxP - 1) / 2; arg1++) {
            double norm = 0;
            double sum_1 =  0;
            double sum_np=0.;
            double sum_nm=0;
            double d1=0;
            double d2=0;
            double d11=0;
            double d22=0;
            double d12=0;
            for (int n2 = -MCp.nMAX; n2 < (MCp.nMAX+1); n2++) {
                for (int n1 = -MCp.nMAX; n1 < (MCp.nMAX+1); n1++) {
                    //std::cout<< "u1: "<< u1<< " u2: "<< u2<< " n1: "<< n1 << " n2: "<< n2<< " norm: "<< norm << std::endl;
                    double u1=dp*arg1 - 2*M_PI*n1;
                    double u2=dp*arg2 - 2*M_PI*n2;
                    //boltz_H = 0.5*(Hp.rho * ((u1*u1) + (u2*u2)) - Hp.nu*(u1-u2)*(u1-u2));
                    double boltz_H = 0.5*(Hp.rho * ((u1*u1) + (u2*u2)) - 2*Hp.nu*(u1*u2));
                    double boltz = exp(-my_beta*boltz_H);
                    norm += boltz;
                    sum_1 += boltz_H*boltz;

                    sum_np+=exp(-beta_np*boltz_H);
                    sum_nm+=exp(-beta_nm*boltz_H);

                    double j1=Hp.rho*u1 - Hp.nu*u2;
                    double j2=Hp.rho*u2 - Hp.nu*u1;

                    d1+= j1*boltz;
                    d2+= j2*boltz;
                    d11+= (Hp.rho - my_beta*j1*j1)*boltz;
                    d22+= (Hp.rho - my_beta*j2*j2)*boltz;
                    d12+= (-Hp.nu - my_beta*j1*j2)*boltz;

                }
            }
            /*Table Villain potential*/
            vil.potential.at(OFFSET_POT + arg1 + MaxP*arg2)= -log(norm)/my_beta;
            fprintf(fVPotential, "%d %d %19.12e \n", arg1, arg2, my_beta*vil.potential[OFFSET_POT + arg1 + MaxP*arg2]);
            /*Calculating beta-derivative of the Villain potential needed for calculating the internal energy (which differs from the energy mean due to the temperature dependence of the Villain model).*/
            vil.upotential.at(OFFSET_POT + (arg1 + MaxP*arg2)) = sum_1/norm;
            vil.potential_bplus.at(OFFSET_POT + arg1 + MaxP*arg2)= -log(sum_np)/beta_np;
            vil.potential_bminus.at(OFFSET_POT + arg1 + MaxP*arg2)= -log(sum_nm)/beta_nm;
            /*Table for the helicity modulus*/
            vil.d1_potential.at(OFFSET_POT + arg1 + MaxP*arg2)+= d1/norm;
            vil.d2_potential.at(OFFSET_POT + arg1 + MaxP*arg2)+= d2/norm;
            vil.d11_potential.at(OFFSET_POT + arg1 + MaxP*arg2)+=d11/norm + my_beta*(d1/norm)*(d1/norm);
            vil.d22_potential.at(OFFSET_POT + arg1 + MaxP*arg2)+=d22/norm + my_beta*(d2/norm)*(d2/norm);
            vil.d12_potential.at(OFFSET_POT + arg1 + MaxP*arg2)+= d12/norm + my_beta*(d1*d2)/(norm*norm);


        }
    }
    fclose(fVPotential);
}




