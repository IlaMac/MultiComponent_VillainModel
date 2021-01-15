//
// Created by ilaria on 2020-12-15.
//

#ifndef VILLAIN_2COMPONENT_PT_H
#define VILLAIN_2COMPONENT_PT_H

#include "main.h"
#include <fstream>
#include "villain_potential.h"

struct PT_parameters{
    /*Parallel Tempering parameters*/
    int np;
    int rank;
    int root=0;
};

struct PTroot_parameters{
    /*Arrays root needs to handle the swaps*/
    std::vector <double> beta;
    std::vector <double> beta_p;
    std::vector <double> beta_m;

    //std::vector <double> All_Energies;
    std::vector <double> E_rank_beta; /* H(x_rank, beta_rank)*/
    std::vector <double> E_rank_betap; /* H(x_rank, beta_{rank +1} )*/
    std::vector <double> E_rank_betam; /* H(x_rank, beta_{rank-1})*/

    std::vector <struct Villain> Villain_beta;

    std::vector <int> ind_to_rank;
    std::vector <int> rank_to_ind;
};

void initialize_PTarrays(struct PT_parameters &PTp, struct PTroot_parameters &PTroot, struct H_parameters &Hp);
void parallel_temp(double &my_E , double &E_betanp, double &E_betanm, double &beta_p, double &beta_m, double &my_beta,  int &my_ind, struct Villain &vil, struct PT_parameters &PTp, struct PTroot_parameters &PTroot);

#endif //VILLAIN_2COMPONENT_PT_H
