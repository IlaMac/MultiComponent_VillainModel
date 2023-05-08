//
// Created by ilaria on 2020-12-17.
//

#ifndef VILLAIN_2COMPONENT_VILLAIN_MC_H
#define VILLAIN_2COMPONENT_VILLAIN_MC_H

#include "constants.h"
#include "main.h"
#include "initialization.h"
#include "villain_potential.h"
#include "rng.h"
#include "class_tic_toc.h"
#include <iostream>
#include <cstring>
#include "o2.h"
#include "measures.h"

void wolff(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta, struct Villain &vil);
void growCluster(int i, int* clusterSpin, const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta, struct Villain &vil);
void wolff_general(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta, struct Villain &vil);
void growCluster_general(int i, int  alpha, int rand_phase, int* clusterSpin, const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta, struct Villain &vil);
void metropolis_villain(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta, struct Villain &vil);
double arg_phase(double x);
int arg(int x, int Max);
int int_arg_phase(int x,  int Max);

#endif //VILLAIN_2COMPONENT_VILLAIN_MC_H
