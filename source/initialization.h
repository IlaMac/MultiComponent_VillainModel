//
// Created by ilaria on 2019-11-13.
//
#pragma once

#include "constants.h"
#include "robust_filesystem.h"
#include "rng.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>

#define Annealing (0) //To be implemented

//static constexpr double C_TWO_PI = 6.2831853071795864769252867665590058;


struct Node{
    /*three spatial dimensions*/
    mutable std::array<double,DIM> A{};
    /*two SC components*/
    mutable std::array<int, NC> Psi{}; /*discretized values of the phase. To be multiplied by 2*M_PI/MaxP*/
};

struct H_parameters{
    /*These values are externally given by an input file*/
    double rho;
    double alpha; //anisotropy factor between the bare stiffness of the two components
    double eta1;
    double eta2;
    double e;
    double h;
    double nu;
    double b_low; //lowest beta for the parallel tempering procedure
    double b_high; // highest beta for the parallel tempering procedure
    int init;
};


struct MC_parameters{
    /*These values are externally given by an input file*/
    int tau; // estimate of the auto-correlation time
    int nmisu; //total number of independent measures
    int n_autosave; //frequency at which intermediate configuration are saved
    //double lbox_theta; //length of the box for the uniform distribution of theta (polar transformation of Psi --> phase)
    double lbox_A; //length of the box for the uniform distribution of dA (transformation of the vector potential)
    int lbox; //length of the box of integer for the extraction of the new single phases
    int lbox_coupled; //length of the box of integer for the extraction of the new coupled phases
    int nMAX; // the sum over integer number in the villain approximation will run from -INT_NMAX to +INT_MAX
};

void initialize_lattice(const std::vector<Node> &Site, const fs::path & directory_read, int RESTART, struct H_parameters &Hp);
void initialize_Hparameters(struct H_parameters &Hp, const fs::path & directory_parameters);
void initialize_MCparameters(struct MC_parameters &MCp, const fs::path & directory_parameters);

