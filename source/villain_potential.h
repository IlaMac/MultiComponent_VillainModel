//
// Created by ilaria on 2020-12-16.
//

#ifndef VILLAIN_2COMPONENT_VILLAIN_POTENTIAL_H
#define VILLAIN_2COMPONENT_VILLAIN_POTENTIAL_H


#include "initialization.h"
#include <fstream>
#include "constants.h"


struct Villain{
    std::vector<double> potential; /*Table of the Villain potential multiplied by beta for each value of \Delta_mu \theta_1(r),  \Delta_mu \theta_2(r) */
    std::vector<double> potential_bplus; /*Table of the Villain potential multiplied by the right neighbouring beta for each value of \Delta_mu \theta_1(r),  \Delta_mu \theta_2(r) */
    std::vector<double> potential_bminus; /*Table of the Villain potential multiplied by the right neighbouring beta for each value of \Delta_mu \theta_1(r),  \Delta_mu \theta_2(r) */
    std::vector<double> upotential; /*Useful to compute efficiently the internal energy*/
    std::vector<double> d1_potential; /*Useful to compute efficiently the helicity modulus*/
    std::vector<double> d2_potential; /*Useful to compute efficiently the helicity modulus*/
    std::vector<double> d11_potential; /*Useful to compute efficiently the helicity modulus*/
    std::vector<double> d22_potential; /*Useful to compute efficiently the helicity modulus*/
    std::vector<double> d12_potential; /*Useful to compute efficiently the helicity modulus*/
};

void init_villain_potentials(double my_beta, double beta_np, double beta_nm, struct Villain &vil, struct H_parameters &Hp, struct MC_parameters &MCp, const fs::path & directory_write);
#endif //VILLAIN_2COMPONENT_VILLAIN_POTENTIAL_H
