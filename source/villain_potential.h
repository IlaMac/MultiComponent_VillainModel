//
// Created by ilaria on 2020-12-16.
//

#ifndef VILLAIN_2COMPONENT_VILLAIN_POTENTIAL_H
#define VILLAIN_2COMPONENT_VILLAIN_POTENTIAL_H


#include "main.h"
#include "initialization.h"
#include <fstream>

#define MaxP (129)

struct Villain{
    double potential[MaxP*MaxP]={0}; /*Table of the Villain potential multiplied by beta for each value of \Delta_mu \theta_1(r),  \Delta_mu \theta_2(r) */
    double potential_bplus[MaxP*MaxP]={0}; /*Table of the Villain potential multiplied by the right neighbouring beta for each value of \Delta_mu \theta_1(r),  \Delta_mu \theta_2(r) */
    double potential_bminus[MaxP*MaxP]={0}; /*Table of the Villain potential multiplied by the right neighbouring beta for each value of \Delta_mu \theta_1(r),  \Delta_mu \theta_2(r) */
    double upotential[MaxP*MaxP]={0}; /*Useful to compute efficiently the internal energy*/
    double d1_potential[MaxP*MaxP]={0}; /*Useful to compute efficiently the helicity modulus*/
    double d2_potential[MaxP*MaxP]={0}; /*Useful to compute efficiently the helicity modulus*/
    double d11_potential[MaxP*MaxP]={0}; /*Useful to compute efficiently the helicity modulus*/
    double d22_potential[MaxP*MaxP]={0}; /*Useful to compute efficiently the helicity modulus*/
    double d12_potential[MaxP*MaxP]={0}; /*Useful to compute efficiently the helicity modulus*/
};

void init_villain_potentials(double my_beta, struct Villain &vil, struct H_parameters &Hp, struct MC_parameters &MCp, const fs::path & directory_write);
void init_villainpotential_nnbeta(double beta_np, double beta_nm, struct Villain &vil,  struct H_parameters &Hp, struct MC_parameters &MCp, const fs::path & directory_write);
#endif //VILLAIN_2COMPONENT_VILLAIN_POTENTIAL_H
