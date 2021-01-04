//
// Created by ilaria on 2019-11-16.
//

#include "measures.h"

void u_internal_energy(struct Measures &mis, struct Villain &vil, struct Node* Site) {

    unsigned int vec;
    int arg_1, arg_2, start=0.5*(MaxP*MaxP-1);;
    double dp=C_TWO_PI/MaxP;
    unsigned int i, ix, iy, iz, nn_i;
    double inv_V=1./N;

    for(iz=0;iz<Lz;iz++){
        for(iy=0;iy<Ly;iy++){
            for(ix=0; ix<Lx; ix++) {
                i = ix + Lx * (iy + iz * Ly);

                for(vec=0; vec<3; vec++) {
                    if (vec == 0) {
                        nn_i = mod(ix + 1, Lx) + Lx * (iy + iz * Ly);
                    }
                    if (vec == 1) {
                        nn_i = ix + Lx * (mod(iy + 1, Ly) + iz * Ly);
                    }
                    if (vec == 2) {
                        nn_i = ix + Lx * (iy + mod(iz + 1, Lz) * Ly);
                    }

                    arg_1 = ARG( (Site[nn_i].Psi[0] - Site[i].Psi[0]), MaxP);
                            //arg((Site[nn_i].Psi[0] - Site[i].Psi[0]));
                    arg_2 = ARG( (Site[nn_i].Psi[1] - Site[i].Psi[1]), MaxP);
                            //arg((Site[nn_i].Psi[1] - Site[i].Psi[1]));
                    mis.U+= vil.upotential[start + arg_1 +MaxP*arg_2];
                }

            }
        }
    }


}

void energy_nn(struct Villain &vil, double &E_betanp, double &E_betanm, struct Node* Site){

    unsigned int vec;
    int arg_1, arg_2, start=0.5*(MaxP*MaxP-1);;
    unsigned int i, ix, iy, iz, nn_i;
    double E_betanp_temp=0., E_betanm_temp=0.;

    for(iz=0;iz<Lz;iz++){
        for(iy=0;iy<Ly;iy++){
            for(ix=0; ix<Lx; ix++){
                i = ix + Lx * (iy + iz * Ly);
                for(vec=0; vec<3; vec++) {
                    if(vec==0){
                        nn_i= mod(ix+1,Lx) + Lx * (iy + iz * Ly);
                    }
                    if(vec==1){
                        nn_i= ix + Lx * (mod(iy+1,Ly) + iz * Ly);
                    }
                    if(vec==2){
                        nn_i= ix + Lx * (iy + mod(iz+1,Lz) * Ly);
                    }

                    arg_1 = ARG( (Site[nn_i].Psi[0] - Site[i].Psi[0]), MaxP);
                            //arg((Site[nn_i].Psi[0] - Site[i].Psi[0]));
                    arg_2 = ARG( (Site[nn_i].Psi[1] - Site[i].Psi[1]), MaxP);
                            //arg((Site[nn_i].Psi[1]- Site[i].Psi[1]));
                    E_betanm_temp+= vil.potential_bminus[start + arg_1 +MaxP*arg_2];
                    E_betanp_temp+= vil.potential_bplus[start + arg_1 +MaxP*arg_2];
                }
            }
        }
    }
    E_betanp=E_betanp_temp;
    E_betanm=E_betanp_temp;
}


void energy(struct Measures &mis, struct Villain &vil, struct Node* Site, double my_beta){

    unsigned int vec;
    int arg_1, arg_2, start=0.5*(MaxP*MaxP-1);;
    unsigned int i, ix, iy, iz, nn_i;

    for(iz=0;iz<Lz;iz++){
        for(iy=0;iy<Ly;iy++){
            for(ix=0; ix<Lx; ix++){
                i = ix + Lx * (iy + iz * Ly);
                for(vec=0; vec<3; vec++) {
                    if(vec==0){
                        nn_i= mod(ix+1,Lx) + Lx * (iy + iz * Ly);
                    }
                    if(vec==1){
                        nn_i= ix + Lx * (mod(iy+1,Ly) + iz * Ly);
                    }
                    if(vec==2){
                        nn_i= ix + Lx * (iy + mod(iz+1,Lz) * Ly);
                    }

                    arg_1 = ARG( (Site[nn_i].Psi[0] - Site[i].Psi[0]), MaxP);
                            //arg((Site[nn_i].Psi[0] - Site[i].Psi[0]));
                    arg_2 = ARG( (Site[nn_i].Psi[1] - Site[i].Psi[1]), MaxP);
                            //arg((Site[nn_i].Psi[1]- Site[i].Psi[1]));
                    mis.E+= vil.potential[start + arg_1 +MaxP*arg_2]/my_beta;

                }
            }
        }
    }

}


void helicity_modulus(double my_beta, struct Measures &mis, struct Villain &vil, struct Node* Site){

    unsigned int vec=0; //I compute the helicity modulus only along one direction x
    int arg_1, arg_2, start=0.5*(MaxP*MaxP-1);;
    double d1=0., d2=0., d11=0., d22=0., d12=0.;
    unsigned int i, ix, iy, iz, nn_i;

    for(iz=0;iz<Lz;iz++){
        for(iy=0;iy<Ly;iy++){
            for(ix=0; ix<Lx; ix++) {
                i = ix + Lx * (iy + iz * Ly);
                nn_i = mod(ix + 1, Lx) + Lx * (iy + iz * Ly);

                arg_1 = ARG( (Site[nn_i].Psi[0] - Site[i].Psi[0]), MaxP);
                        //arg((Site[nn_i].Psi[0] - Site[i].Psi[0]));
                arg_2 = ARG( (Site[nn_i].Psi[1] - Site[i].Psi[1]), MaxP);
                        //arg((Site[nn_i].Psi[1] - Site[i].Psi[1]));

                d1+= vil.d1_potential[start + arg_1 +MaxP*arg_2];
                d2+= vil.d2_potential[start + arg_1 +MaxP*arg_2];
                d11+= vil.d11_potential[start + arg_1 +MaxP*arg_2];
                d22+= vil.d22_potential[start + arg_1 +MaxP*arg_2];
                d12+= vil.d12_potential[start + arg_1 +MaxP*arg_2];

            }
        }
    }

    mis.DH_Ddi[0]=d1;
    mis.DH_Ddi[1]=d2;
    mis.D2H_Dd2i[0]=d11;
    mis.D2H_Dd2i[1]=d22;
    mis.D2H_Dd12=d12;

}

//TO CHANGE!!!!!!!
//void dual_stiffness(struct Measures &mis, struct H_parameters &Hp, struct Node* Site){
//
//    double qx_min=C_TWO_PI/(Lx);
//    double invNorm= 1./((C_TWO_PI)*(C_TWO_PI)*N);
//    unsigned int i, ix, iy, iz;
//    double Re_rhoz=0.;
//    double Im_rhoz=0.;
//    double Dx_Ay, Dy_Ax;
//
//    for(ix=0; ix<Lx;ix++){
//        for(iy=0; iy<Ly;iy++){
//            for(iz=0; iz<Lz;iz++){
//                i=ix +Lx*(iy+Ly*iz);
//                Dx_Ay=(Site[mod(ix + 1, Lx) +Lx*(iy+Ly*iz)].A[1]- Site[i].A[1])/Hp.h;
//                Dy_Ax=(Site[ix +Lx*(mod(iy + 1, Ly) +Ly*iz) ].A[0]- Site[i].A[0])/Hp.h;
//
//                Re_rhoz+=(cos((double)qx_min*ix)*(Dx_Ay -Dy_Ax));
//                Im_rhoz+=(sin((double)qx_min*ix)*(Dx_Ay -Dy_Ax));
//            }
//        }
//    }
//    mis.d_rhoz=invNorm*((Re_rhoz*Re_rhoz) +(Im_rhoz*Im_rhoz));
//}

//DESIGNED FOR 3 COMPONENTs
void magnetization(struct Measures &mis, struct Node* Site){
    //The Ising parameter m(x,y)=+/-1 indicates the chirality of the three phases. If the phases are ordered as: phi_1, phi_2, phi_3 then m=1; otherwise if the order is phi_1, phi_3, phi_2 then m=-1.
    unsigned ix, iy, iz, i;

    double phi_shifted_1=0.;
    double phi_shifted_2=0.;

    for(iz=0; iz<Lz;iz++) {
        for (iy = 0; iy < Ly; iy++) {
            for (ix = 0; ix < Lx; ix++) {
                i=ix +Lx*(iy+Ly*iz);

                phi_shifted_1= Site[i].Psi[1] - Site[i].Psi[0];
                while(phi_shifted_1 >= C_TWO_PI){
                        phi_shifted_1-= C_TWO_PI;}
                while(phi_shifted_1< 0){
                        phi_shifted_1+=C_TWO_PI;}

                phi_shifted_2= Site[i].Psi[2] - Site[i].Psi[0];
                while(phi_shifted_2 >= C_TWO_PI){
                        phi_shifted_2-= C_TWO_PI;}
                while(phi_shifted_2< 0){
                        phi_shifted_2+=C_TWO_PI;}

                if(phi_shifted_1>phi_shifted_2){
                    mis.m+=1;
                        }else if(phi_shifted_1<phi_shifted_2){
                    mis.m+=(-1);
                }
            }
        }
    }
    mis.m=mis.m/N;
}

void magnetization_singlephase(struct Measures &mis, struct Node* Site, double my_beta){
    //The Ising parameter m(x,y)=+/-1 indicates the chirality of the three phases. If the phases are ordered as: phi_1, phi_2, phi_3 then m=1; otherwise if the order is phi_1, phi_3, phi_2 then m=-1.
    unsigned ix, iy, iz, i, alpha;
    double dp=C_TWO_PI/MaxP;
    double cos_phi[NC]={0}, sin_phi[NC]={0};
    double inv_N=1./N;
    for(iz=0; iz<Lz;iz++) {
        for (iy = 0; iy < Ly; iy++) {
            for (ix = 0; ix < Lx; ix++) {
                i=ix +Lx*(iy+Ly*iz);
                for(alpha=0; alpha<NC; alpha++){
                    cos_phi[alpha]+= cos(Site[i].Psi[alpha]*dp);
                    sin_phi[alpha]+= sin(Site[i].Psi[alpha]*dp);
                    //std::cout << "alpha: "<< alpha <<  " phase: "<< Site[i].Psi[alpha]*dp << "int phase: " << Site[i].Psi[alpha]<< std::endl;
                }
            }
        }
    }

    for(alpha=0; alpha<NC; alpha++) {
        cos_phi[alpha]*=inv_N;
        sin_phi[alpha]*=inv_N;
        mis.m_phase[alpha] = (cos_phi[alpha]*cos_phi[alpha]) + (sin_phi[alpha]*sin_phi[alpha]);
    }

 //std::cout << "my beta"<< my_beta<<  " mphase1: "<< mis.m_phase[0] << " mphase2: " << mis.m_phase[1]<< std::endl;


}

void save_lattice(struct Node* Site, const fs::path & directory_write, std::string configuration){

    std::string sPsi;
    std::string sA;
    sPsi= std::string("Psi_")+ configuration + std::string(".bin");
    sA= std::string("A_")+ configuration + std::string(".bin");
    fs::path psi_init_file = directory_write / sPsi;
    fs::path a_init_file = directory_write / sA;

    FILE *fPsi= nullptr;
    FILE *fA= nullptr;
    unsigned int i=0;

    if((fPsi=fopen(psi_init_file.c_str(), "w"))) {
        for (i = 0; i < N; i++) {
            fwrite(Site[i].Psi, sizeof(int), NC, fPsi);
        }
        fclose(fPsi);
    }

    if((fA=fopen(a_init_file.c_str(), "w"))) {
    	for (i = 0; i < N; i++) {
            fwrite(Site[i].A, sizeof(double), 3, fA);
	}
        fclose(fA);
    }

}

