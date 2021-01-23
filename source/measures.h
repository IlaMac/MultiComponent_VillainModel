//
// Created by ilaria on 2019-11-16.
//

#ifndef MEASURES_H
#define MEASURES_H

#include "constants.h"
#include "main.h"
#include <fstream>
#include "villain_potential.h"
#include "villain_MC.h"

struct Measures{

/*******THE ORDER SHOULD BE THE SAME OF THE H5T INSERT REGISTRATION**************/

    double E=0.; //Energy
    double U=0.;
    double m=0.; //magnetization (for the phase chirality of the three components
    //Binder cumulant U=<m⁴>/(3*<m²>²)
    double m_phase[NC]={0}; //magnetization of the single component phase
    double d_rhoz=0.; //Dual stiffness along z

    double DH_Ddi[NC]={0}; //1st derivative in the twisted phase of the i component
    double D2H_Dd2i[NC]={0}; //2nd derivative in the twisted phase of the i component
    double D2H_Dd12=0.; //2nd mixed derivative in the twisted phases of the component 1 and 2

    int my_rank = 0;
    void reset(){
        *this = Measures();
    }
};


void energy(struct Measures &mis, struct Villain &vil, const std::vector<Node> &Site, double my_beta);
void energy_nn(struct Villain &vil, double &E_betanp, double &E_betanm, const std::vector<Node> &Site);
void u_internal_energy(struct Measures &mis, struct Villain &vil, const std::vector<Node> &Site);
void helicity_modulus(struct Measures &mis, struct Villain &vil, const std::vector<Node> &Site);
void magnetization(struct Measures &mis, const std::vector<Node> &Site);
void magnetization_singlephase(struct Measures &mis, const std::vector<Node> &Site, double my_beta);
void save_lattice(const std::vector<Node> &Site, const fs::path & directory_write, std::string configuration);
//void dual_stiffness(struct Measures &mis, struct H_parameters &Hp, const std::vector<Node> &Site);

#endif //MEASURES_H
