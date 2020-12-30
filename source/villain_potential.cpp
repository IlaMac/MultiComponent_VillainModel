//
// Created by ilaria on 2020-12-16.
//

#include "villain_potential.h"

void init_villain_potentials(double my_beta, struct Villain &vil,  struct H_parameters &Hp, struct MC_parameters &MCp) {

    int n1, n2, arg1, arg2, start=0.5*(MaxP*MaxP-1);
    double u1, u2, sum_1, sum, norm, boltz, boltz_H;
    double j1, j2;
    double d1, d2, d11, d12, d22;
    double dp=C_TWO_PI/MaxP;
    double inv_beta=1./my_beta;

    for (arg2 = -(MaxP - 1) / 2; arg2 <= (MaxP - 1) / 2; arg2++) {
        for (arg1 = -(MaxP - 1) / 2; arg1 <= (MaxP - 1) / 2; arg1++) {
            sum = 0;
            norm = 0;
            sum_1 =  0;
            d1=0;
            d2=0;
            d11=0;
            d22=0;
            d12=0;
            for (n2 = -MCp.nMAX; n2 < (MCp.nMAX+1); n2++) {
                for (n1 = -MCp.nMAX; n1 < (MCp.nMAX+1); n1++) {
                    u1=dp*arg1 - C_TWO_PI*n1;
                    u2=dp*arg2 - C_TWO_PI*n2;

                    sum+=exp(-0.5*my_beta*(Hp.rho * (u1*u1 +u2*u2) + Hp.nu*(u1*u2) ));

                    boltz_H = 0.5*(Hp.rho * (u1*u1 +u2*u2) + Hp.nu*(u1*u2) );
                    boltz = exp(-my_beta*boltz_H);
                    norm += boltz;
                    sum_1 += boltz_H*boltz;

                    j1=Hp.rho*u1 + 0.5*Hp.nu*u2;
                    j2=Hp.rho*u2 + 0.5*Hp.nu*u1;

                    d1+= j1*boltz;
                    d2+= j2*boltz;
                    d11+= (Hp.rho - my_beta*j1*j1)*boltz;
                    d22+= (Hp.rho - my_beta*j2*j2)*boltz;
                    d12+= (0.5*Hp.nu - my_beta*j1*j2)*boltz;

                }
            }
            /*Table Villain potential*/
            vil.potential[start + arg1 + MaxP*arg2]= -inv_beta*log(sum);
            /*Table for the helicity modulus*/
            vil.d1_potential[start + arg1 + MaxP*arg2]+= d1/norm;
            vil.d2_potential[start + arg1 + MaxP*arg2]+= d2/norm;
            vil.d11_potential[start + arg1 + MaxP*arg2]+=d11/norm + my_beta*(d1/norm)*(d1/norm);
            vil.d22_potential[start + arg1 + MaxP*arg2]+=d22/norm + my_beta*(d2/norm)*(d2/norm);
            vil.d12_potential[start + arg1 + MaxP*arg2]+= d12/norm + my_beta*(d1*d2)/(norm*norm);
            /*Calculating beta-derivative of the Villain potential needed for calculating the internal energy (which differs from the energy mean due to the temperature dependence of the Villain model).*/
            vil.upotential[start + (arg1 + MaxP*arg2)] = sum_1/norm;
        }
    }

}

void init_villainpotential_nnbeta(double beta_np, double beta_nm, struct Villain &vil,  struct H_parameters &Hp, struct MC_parameters &MCp) {

    int n1, n2, arg1, arg2, start=0.5*(MaxP*MaxP-1);
    double sum_np, sum_nm, u1, u2;
    double dp=C_TWO_PI/MaxP;
    double inv_betanp=1./beta_np;
    double inv_betanm=1./beta_nm;

    for (arg2 = -(MaxP - 1) / 2; arg2 <= (MaxP - 1) / 2; arg2++) {
        for (arg1 = -(MaxP - 1) / 2; arg1 <= (MaxP - 1) / 2; arg1++) {
            sum_np = 0;
            sum_nm = 0;
            for (n2 = -MCp.nMAX; n2 < (MCp.nMAX+1); n2++) {
                for (n1 = -MCp.nMAX; n1 < (MCp.nMAX+1); n1++) {
                    u1=dp*arg1 - C_TWO_PI*n1;
                    u2=dp*arg2 - C_TWO_PI*n2;
                    sum_np+=exp(-0.5*beta_np*(Hp.rho * (u1*u1 +u2*u2) + Hp.nu*(u1*u2) ));
                    sum_nm+=exp(-0.5*beta_nm*(Hp.rho * (u1*u1 +u2*u2) + Hp.nu*(u1*u2) ));
                }
            }
            vil.potential_bplus[start + arg1 + MaxP*arg2]= -inv_betanp*log(sum_np);
            vil.potential_bminus[start + arg1 + MaxP*arg2]= -inv_betanm*log(sum_nm);
        }
    }
}



