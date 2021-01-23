//
// Created by ilaria on 2019-11-16.
//

#include "measures.h"

void u_internal_energy(struct Measures &mis, struct Villain &vil, const std::vector<Node> &Site) {

    unsigned int vec;
    int arg_1, arg_2;
    int i, ix, iy, iz, nn_i;

    for(iz=0;iz<Lz;iz++){
        for(iy=0;iy<Ly;iy++){
            for(ix=0; ix<Lx; ix++) {
                i = ix + Lx * (iy + iz * Ly);

                for(vec=0; vec<3; vec++) {
                    if (vec == 0) {
                        nn_i = mod(ix + 1, (int)Lx) + Lx * (iy + iz * Ly);
                    }
                    if (vec == 1) {
                        nn_i = ix + Lx * (mod(iy + 1, (int)Ly) + iz * Ly);
                    }
                    if (vec == 2) {
                        nn_i = ix + Lx * (iy + mod(iz + 1, (int)Lz) * Ly);
                    }

                    arg_1 = arg((Site.at(nn_i).Psi[0] - Site.at(i).Psi[0]), MaxP);
                    arg_2 = arg((Site[nn_i].Psi[1] - Site.at(i).Psi[1]), MaxP);
                    mis.U+= vil.upotential.at(OFFSET_POT + arg_1 +MaxP*arg_2);
                }

            }
        }
    }


}

void energy_nn(struct Villain &vil, double &E_betanp, double &E_betanm, const std::vector<Node> &Site){

    unsigned int vec;
    int arg_1, arg_2;
    int i, ix, iy, iz, nn_i;
    double E_betanp_temp=0., E_betanm_temp=0.;

    for(iz=0;iz<Lz;iz++){
        for(iy=0;iy<Ly;iy++){
            for(ix=0; ix<Lx; ix++){
                i = ix + Lx * (iy + iz * Ly);
                for(vec=0; vec<3; vec++) {
                    if(vec==0){
                        nn_i= mod(ix+1,(int)Lx) + Lx * (iy + iz * Ly);
                    }
                    if(vec==1){
                        nn_i= ix + Lx * (mod(iy+1,(int)Ly) + iz * Ly);
                    }
                    if(vec==2){
                        nn_i= ix + Lx * (iy + mod(iz+1,(int)Lz) * Ly);
                    }
                    arg_1 = arg((Site[nn_i].Psi[0] - Site.at(i).Psi[0]), MaxP);
                    arg_2 = arg((Site[nn_i].Psi[1]- Site.at(i).Psi[1]), MaxP);
                    E_betanm_temp+= vil.potential_bminus.at(OFFSET_POT + arg_1 +MaxP*arg_2);
                    E_betanp_temp+= vil.potential_bplus.at(OFFSET_POT + arg_1 +MaxP*arg_2);
                }
            }
        }
    }
    E_betanp=E_betanp_temp;
    E_betanm=E_betanp_temp;
}


void energy(struct Measures &mis, struct Villain &vil, const std::vector<Node> &Site){

    unsigned int vec;
    int arg_1, arg_2;
    int i, ix, iy, iz, nn_i;

    for(iz=0;iz<Lz;iz++){
        for(iy=0;iy<Ly;iy++){
            for(ix=0; ix<Lx; ix++){
                i = ix + Lx * (iy + iz * Ly);
                for(vec=0; vec<3; vec++) {
                    if(vec==0){
                        nn_i= mod(ix+1, (int)Lx) + Lx * (iy + iz * Ly);
                    }
                    if(vec==1){
                        nn_i= ix + Lx * (mod(iy+1,(int)Ly) + iz * Ly);
                    }
                    if(vec==2){
                        nn_i= ix + Lx * (iy + mod(iz+1,(int)Lz) * Ly);
                    }

                    arg_1 = arg((Site[nn_i].Psi[0] - Site.at(i).Psi[0]), MaxP);
                    arg_2 = arg((Site[nn_i].Psi[1]- Site.at(i).Psi[1]), MaxP);
                    mis.E+= vil.potential.at(OFFSET_POT + arg_1 +MaxP*arg_2);

                }
            }
        }
    }

}


void helicity_modulus(struct Measures &mis, struct Villain &vil, const std::vector<Node> &Site){

    //I compute the helicity modulus only along one direction x
    int arg_1, arg_2;
    double d1=0., d2=0., d11=0., d22=0., d12=0.;
    int i, ix, iy, iz, nn_i;

    for(iz=0;iz<Lz;iz++){
        for(iy=0;iy<Ly;iy++){
            for(ix=0; ix<Lx; ix++) {
                i = ix + Lx * (iy + iz * Ly);
                nn_i = mod(ix + 1, Lx) + Lx * (iy + iz * Ly);

                arg_1 = arg((Site[nn_i].Psi[0] - Site.at(i).Psi[0]), MaxP);
                arg_2 = arg((Site[nn_i].Psi[1] - Site.at(i).Psi[1]), MaxP);

                d1+= vil.d1_potential.at(OFFSET_POT + arg_1 +MaxP*arg_2);
                d2+= vil.d2_potential.at(OFFSET_POT + arg_1 +MaxP*arg_2);
                d11+= vil.d11_potential.at(OFFSET_POT + arg_1 +MaxP*arg_2);
                d22+= vil.d22_potential.at(OFFSET_POT + arg_1 +MaxP*arg_2);
                d12+= vil.d12_potential.at(OFFSET_POT + arg_1 +MaxP*arg_2);
            }
        }
    }

    mis.DH_Ddi[0]=d1;
    mis.DH_Ddi[1]=d2;
    mis.D2H_Dd2i[0]=d11;
    mis.D2H_Dd2i[1]=d22;
    mis.D2H_Dd12=d12;

}

////DESIGNED FOR 3 COMPONENTs
//void magnetization(struct Measures &mis, const std::vector<Node> &Site){
//    //The Ising parameter m(x,y)=+/-1 indicates the chirality of the three phases. If the phases are ordered as: phi_1, phi_2, phi_3 then m=1; otherwise if the order is phi_1, phi_3, phi_2 then m=-1.
//    unsigned ix, iy, iz, i;
//
//    double phi_shifted_1=0.;
//    double phi_shifted_2=0.;
//
//    for(iz=0; iz<Lz;iz++) {
//        for (iy = 0; iy < Ly; iy++) {
//            for (ix = 0; ix < Lx; ix++) {
//                i=ix +Lx*(iy+Ly*iz);
//
//                phi_shifted_1= Site.at(i).Psi[1] - Site.at(i).Psi[0];
//                while(phi_shifted_1 >= 2*M_PI){
//                        phi_shifted_1-= 2*M_PI;}
//                while(phi_shifted_1< 0){
//                        phi_shifted_1+=2*M_PI;}
//
//                phi_shifted_2= Site.at(i).Psi[2] - Site.at(i).Psi[0];
//                while(phi_shifted_2 >= 2*M_PI){
//                        phi_shifted_2-= 2*M_PI;}
//                while(phi_shifted_2< 0){
//                        phi_shifted_2+=2*M_PI;}
//
//                if(phi_shifted_1>phi_shifted_2){
//                    mis.m+=1;
//                        }else if(phi_shifted_1<phi_shifted_2){
//                    mis.m+=(-1);
//                }
//            }
//        }
//    }
//    mis.m=mis.m/N;
//}

void magnetization_singlephase(struct Measures &mis, const std::vector<Node> &Site){
    //The Ising parameter m(x,y)=+/-1 indicates the chirality of the three phases. If the phases are ordered as: phi_1, phi_2, phi_3 then m=1; otherwise if the order is phi_1, phi_3, phi_2 then m=-1.
    unsigned alpha;
    double cos_phi[NC]={0}, sin_phi[NC]={0};
    double inv_N=1./N;

    for(auto & s : Site) {
        for (alpha = 0; alpha < NC; alpha++) {
            cos_phi[alpha] += cos(s.Psi[alpha] * dp);
            sin_phi[alpha] += sin(s.Psi[alpha] * dp);
        }
    }

    for(alpha=0; alpha<NC; alpha++) {
        cos_phi[alpha]*=inv_N;
        sin_phi[alpha]*=inv_N;
        mis.m_phase[alpha] = (cos_phi[alpha]*cos_phi[alpha]) + (sin_phi[alpha]*sin_phi[alpha]);
    }

 //std::cout << "my beta"<< my_beta<<  " mphase1: "<< mis.m_phase[0] << " mphase2: " << mis.m_phase[1]<< std::endl;


}

void save_lattice(const std::vector<Node> &Site, const fs::path & directory_write, const std::string& configuration){

    std::string sPsi;
    std::string sA;
    sPsi= std::string("Psi_")+ configuration + std::string(".bin");
    sA= std::string("A_")+ configuration + std::string(".bin");
    fs::path psi_init_file = directory_write / sPsi;
    fs::path a_init_file = directory_write / sA;

    FILE *fPsi= nullptr;
    FILE *fA= nullptr;

    if((fPsi=fopen(psi_init_file.c_str(), "w"))) {
        for (auto & s: Site) {
            fwrite(s.Psi.data(), sizeof(int), NC, fPsi);
        }
        fclose(fPsi);
    }

    if((fA=fopen(a_init_file.c_str(), "w"))) {
        for (auto & s: Site) {
            fwrite(s.A.data(), sizeof(double), 3, fA);
	}
        fclose(fA);
    }

}

