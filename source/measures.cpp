//
// Created by ilaria on 2019-11-16.
//

#include "measures.h"

void energy(struct Measures &mis, struct H_parameters &Hp, struct MC_parameters &MCp, double my_beta, struct Node* Site){

    unsigned int vec;
    int n_1, n_2;
    unsigned int i, ix, iy, iz, nn_i;
    double u0_1, u0_2, u_1, u_2;
    double local_S=0., sum_n1n2=0., boltz_h=0.;
    double inv_beta=1./my_beta;

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

                    u0_1= Site[nn_i].Psi[0].t - Site[i].Psi[0].t;
                    u0_2= Site[nn_i].Psi[1].t - Site[i].Psi[1].t;

                    for (n_1 = -(MCp.nMAX); n_1 < (MCp.nMAX+1); n_1++) {
                        for (n_2 = -(MCp.nMAX); n_2 < (MCp.nMAX+1); n_2++) {

                            u_1= u0_1 - C_TWO_PI*n_1;
                            u_2= u0_2 - C_TWO_PI*n_2;

                            local_S = 0.5 * my_beta * (Hp.rho * (u_1*u_1 + u_2*u_2) + Hp.nu*(u_1*u_2) );

                            sum_n1n2+= exp(-local_S);
                        }
                    }

                    boltz_h-=inv_beta*log(sum_n1n2);
                }
            }
        }
    }

    mis.E=boltz_h; // As soon as I will include the lattice spacing h --> h3*boltz_h
}


void helicity_modulus(struct Measures &mis, struct H_parameters &Hp, struct MC_parameters &MCp, double my_beta, struct Node* Site){

    unsigned int vec=0; //I compute the helicity modulus only along one direction x
    int n_1, n_2;
    unsigned int i, ix, iy, iz, nn_i;
    double u_1, u_2, u0_1, u0_2;
    double j1, j2;
    double d1=0., d2=0., d11=0., d12=0., d22=0.;
    double norm=0., boltz=0.;
    double inv_beta=1./my_beta;
    double inv_V=1./N;

    for(iz=0;iz<Lz;iz++){
        for(iy=0;iy<Ly;iy++){
            for(ix=0; ix<Lx; ix++) {
                i = ix + Lx * (iy + iz * Ly);
                nn_i = mod(ix + 1, Lx) + Lx * (iy + iz * Ly);

                u0_1= Site[nn_i].Psi[0].t - Site[i].Psi[0].t;
                u0_2= Site[nn_i].Psi[1].t - Site[i].Psi[1].t;

                for (n_1 = -MCp.nMAX; n_1 < (MCp.nMAX +1); n_1++) {
                    for (n_2 = -MCp.nMAX; n_2 < (MCp.nMAX+1); n_2++) {
                        u_1 = u0_1 - C_TWO_PI * n_1;
                        u_2 = u0_2 - C_TWO_PI * n_2;

                        j1=Hp.rho*u_1 + Hp.nu*u_2;
                        j2=Hp.rho*u_2 + Hp.nu*u_1;
                        boltz=exp(-(0.5 * my_beta * (Hp.rho * (u_1*u_1 + u_2*u_2) + Hp.nu*(u_1*u_2) )));

                        norm+=boltz;
                        d1+= j1*boltz;
                        d2+= j2*boltz;
                        d11+= (Hp.rho - my_beta*j1*j1)*boltz;
                        d22+= (Hp.rho - my_beta*j2*j2)*boltz;
                        d12+= (Hp.nu - my_beta*j1*j2)*boltz;

                    }
                }

                mis.DH_Ddi[0]+= d1/norm;
                mis.DH_Ddi[1]+= d2/norm;
                mis.D2H_Dd2i[0]+=d11/norm + my_beta*(d1/norm)*(d1/norm);
                mis.D2H_Dd2i[1]+=d22/norm + my_beta*(d2/norm)*(d2/norm);
                mis.D2H_Dd12+= d12/norm + my_beta*(d1*d2)/(norm*norm);

            }
        }
    }

    mis.DH_Ddi[0]*=inv_V;
    mis.DH_Ddi[1]*=inv_V;
    mis.D2H_Dd2i[0]*=inv_V;
    mis.D2H_Dd2i[1]*=inv_V;
    mis.D2H_Dd12*=inv_V;

}

//TO CHANGE!!!!!!!
void dual_stiffness(struct Measures &mis, struct H_parameters &Hp, struct Node* Site){

    double qx_min=C_TWO_PI/(Lx);
    double invNorm= 1./((C_TWO_PI)*(C_TWO_PI)*N);
    unsigned int i, ix, iy, iz;
    double Re_rhoz=0.;
    double Im_rhoz=0.;
    double Dx_Ay, Dy_Ax;

    for(ix=0; ix<Lx;ix++){
        for(iy=0; iy<Ly;iy++){
            for(iz=0; iz<Lz;iz++){
                i=ix +Lx*(iy+Ly*iz);
                Dx_Ay=(Site[nn(i, 0, 1)].A[1]- Site[i].A[1])/Hp.h;
                Dy_Ax=(Site[nn(i, 1, 1)].A[0]- Site[i].A[0])/Hp.h;

                Re_rhoz+=(cos((double)qx_min*ix)*(Dx_Ay -Dy_Ax));
                Im_rhoz+=(sin((double)qx_min*ix)*(Dx_Ay -Dy_Ax));
            }
        }
    }
    mis.d_rhoz=invNorm*((Re_rhoz*Re_rhoz) +(Im_rhoz*Im_rhoz));
}

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

                phi_shifted_1= Site[i].Psi[1].t - Site[i].Psi[0].t;
                while(phi_shifted_1 >= C_TWO_PI){
                        phi_shifted_1-= C_TWO_PI;}
                while(phi_shifted_1< 0){
                        phi_shifted_1+=C_TWO_PI;}

                phi_shifted_2= Site[i].Psi[2].t - Site[i].Psi[0].t;
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

void magnetization_singlephase(struct Measures &mis, struct Node* Site){
    //The Ising parameter m(x,y)=+/-1 indicates the chirality of the three phases. If the phases are ordered as: phi_1, phi_2, phi_3 then m=1; otherwise if the order is phi_1, phi_3, phi_2 then m=-1.
    unsigned ix, iy, iz, i, alpha;
    double cos_phi[NC]={0}, sin_phi[NC]={0};
    double inv_N=1./N;
    for(iz=0; iz<Lz;iz++) {
        for (iy = 0; iy < Ly; iy++) {
            for (ix = 0; ix < Lx; ix++) {
                i=ix +Lx*(iy+Ly*iz);
                for(alpha=0; alpha<NC; alpha++){
                    cos_phi[alpha]+= cos(Site[i].Psi[alpha].t);
                    sin_phi[alpha]+= sin(Site[i].Psi[alpha].t);
                }
            }
        }
    }

    for(alpha=0; alpha<NC; alpha++) {
        cos_phi[alpha]*=inv_N;
        sin_phi[alpha]*=inv_N;
        mis.m_phase[alpha] = (cos_phi[alpha]*cos_phi[alpha]) + (sin_phi[alpha]*sin_phi[alpha]);
    }
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
            fwrite(Site[i].Psi, sizeof(struct O2), 3, fPsi);
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

