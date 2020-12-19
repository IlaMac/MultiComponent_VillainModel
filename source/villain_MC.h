//
// Created by ilaria on 2020-12-17.
//

#ifndef VILLAIN_2COMPONENT_VILLAIN_MC_H
#define VILLAIN_2COMPONENT_VILLAIN_MC_H

#include "main.h"
#include "initialization.h"
#include "villain_potential.h"
#include "rng.h"
#include "class_tic_toc.h"
#include <iostream>
#include <cstring>
#include "o2.h"



#undef  ARG
#define ARG(diff,N) ((diff) + (N)*((2*(diff)<(-(N)))*(1-(2*(diff)+(N)+1)/(2*(N))) - (2*(diff)>(N))*(1+(2*(diff)-(N)-1)/(2*(N)))))


void metropolis_villain(struct Node* Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta, struct Villain &vil);
double arg_phase(double x);
int arg(int x);
int int_arg_phase(int x);

#endif //VILLAIN_2COMPONENT_VILLAIN_MC_H
