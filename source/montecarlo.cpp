#include "montecarlo.h"
#include "main.h"
#include "rng.h"
#include "class_tic_toc.h"

void metropolis( struct Node* Site, struct MC_parameters &MCp, struct H_parameters &Hp,  double my_beta){

    double d_theta, d_A, rand;
    unsigned int ix, iy, iz, alpha, vec, i;
    double acc_rate=0.5, acc_A=0., acc_theta=0.; // acc_rho=0.,
    struct O2 NewPsi{};
    struct O2 OldPsi{};
    double NewA, OldA;
    double newE, oldE, minus_deltaE;
    double h3=(Hp.h*Hp.h*Hp.h);
    class_tic_toc t_localHtheta(true,5,"local_Htheta");
    class_tic_toc t_localHA(true,5,"local_HA");


    for (iz= 0; iz < Lz; iz++) {
        for (iy = 0; iy < Ly; iy++) {
            for (ix = 0; ix < Lx; ix++) {
                i = ix + Lx * (iy + iz * Ly);

                /*******PHASE ONLY UPDATE**************/
                for (alpha = 0; alpha < NC; alpha++) {
                    OldPsi = Site[i].Psi[alpha];
                    //t_localHtheta.tic();
                    oldE = HVillan_old(i, alpha, my_beta, Hp, MCp, Site);
                    //t_localHtheta.toc();
                    d_theta = rn::uniform_real_box(-MCp.lbox_theta, MCp.lbox_theta);
                    NewPsi.t = fmod(OldPsi.t + d_theta, C_TWO_PI);
                    NewPsi.r=1;
                    //t_localHtheta.tic();
                    newE = HVillan_new(NewPsi, i, alpha, my_beta, Hp, MCp, Site);
                    //t_localHtheta.toc();
                    minus_deltaE =  (oldE - newE);
                    // h3 * (oldE - newE); when I will include properly the lattice spacing which in this case is taken to be just 1
                    if (minus_deltaE > 0) {
                        Site[i].Psi[alpha] = NewPsi;
                        acc_theta++;
                    } else {
                        rand = rn::uniform_real_box(0, 1);
                        //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
                        if (rand < exp( minus_deltaE)) {
                            Site[i].Psi[alpha] = NewPsi;
                            acc_theta++;
                        }
                    }
                }
            }
        }
    }
    acc_theta=(double) acc_theta/(3*N);
    //acc_A=(double) acc_A/(3*N);
    MCp.lbox_theta= MCp.lbox_theta*(0.5*(acc_theta/acc_rate)+0.5);
    //MCp.lbox_A= MCp.lbox_A*(0.5*(acc_A/acc_rate)+0.5);
}

double HVillan_old(unsigned int position, unsigned int alpha, double my_beta, struct H_parameters &Hp, struct MC_parameters &MCp, struct Node* Site){

    unsigned int vec;
    int n_1, n_2;
    unsigned int ix, iy, iz, nn_ip, nn_im;
    double u0_1p, u0_2p, u0_1m, u0_2m, u_1, u_2;
    double local_S=0., sum_n1n2_plus=0., sum_n1n2_minus=0., boltz_h=0.;
    double inv_beta=1./my_beta;

    ix=position%Lx;
    iy=(position/Lx)%Ly;
    iz=(position/(Lx*Ly));

    for(vec=0; vec<3; vec++) {
        if(vec==0){
            nn_ip= mod(ix+1,Lx) + Lx * (iy + iz * Ly);
            nn_im= mod(ix-1,Lx) + Lx * (iy + iz * Ly);
        }
        if(vec==1){
            nn_ip= ix + Lx * (mod(iy+1,Ly) + iz * Ly);
            nn_im= ix + Lx * (mod(iy-1,Ly) + iz * Ly);
        }
        if(vec==2){
            nn_ip= ix + Lx * (iy + mod(iz+1,Lz) * Ly);
            nn_im= ix + Lx * (iy + mod(iz-1,Lz) * Ly);
        }

        u0_1p= Site[nn_ip].Psi[0].t - Site[position].Psi[0].t;
        u0_2p= Site[nn_ip].Psi[1].t - Site[position].Psi[1].t;

        u0_1m= -Site[nn_im].Psi[0].t + Site[position].Psi[0].t;
        u0_2m= -Site[nn_im].Psi[1].t + Site[position].Psi[1].t;

        for (n_1 = (-MCp.nMAX); n_1 < (MCp.nMAX +1); n_1++) {
            for (n_2 = (-MCp.nMAX); n_2 < (MCp.nMAX +1); n_2++) {

                u_1= u0_1p - C_TWO_PI*n_1;
                u_2= u0_2p - C_TWO_PI*n_2;

                local_S = 0.5 * my_beta * (Hp.rho * (u_1*u_1 + u_2*u_2) + Hp.nu*(u_1*u_2) );
                sum_n1n2_plus+= exp(-local_S);

                u_1= u0_1m - C_TWO_PI*n_1;
                u_2= u0_2m - C_TWO_PI*n_2;

                local_S = 0.5 * my_beta * (Hp.rho * (u_1*u_1 + u_2*u_2) + Hp.nu*(u_1*u_2) );
                sum_n1n2_minus+= exp(-local_S);
            }
        }

        boltz_h-=log(sum_n1n2_minus);
        boltz_h-=log(sum_n1n2_plus);

    }


    return boltz_h;

}

double HVillan_new(struct O2 Psi, unsigned int position, unsigned int alpha,  double my_beta, struct H_parameters &Hp, struct MC_parameters &MCp, struct Node* Site){

    unsigned int vec;
    int n_1, n_2;
    unsigned int ix, iy, iz, nn_ip, nn_im;
    double theta_1=0., theta_2=0.;
    double u0_1p, u0_2p, u0_1m, u0_2m, u_1, u_2;
    double local_S=0., sum_n1n2_plus=0., sum_n1n2_minus=0., boltz_h=0.;
    double inv_beta=1./my_beta;

    ix=position%Lx;
    iy=(position/Lx)%Ly;
    iz=(position/(Lx*Ly));

    if(alpha==0){
        theta_1=Psi.t;
        theta_2=Site[position].Psi[1].t;
    }
    else if(alpha==1){
        theta_1=Site[position].Psi[0].t;
        theta_2=Psi.t;
    }

    for(vec=0; vec<3; vec++) {
        if(vec==0){
            nn_ip= mod(ix+1,Lx) + Lx * (iy + iz * Ly);
            nn_im= mod(ix-1,Lx) + Lx * (iy + iz * Ly);
        }
        if(vec==1){
            nn_ip= ix + Lx * (mod(iy+1,Ly) + iz * Ly);
            nn_im= ix + Lx * (mod(iy-1,Ly) + iz * Ly);
        }
        if(vec==2){
            nn_ip= ix + Lx * (iy + mod(iz+1,Lz) * Ly);
            nn_im= ix + Lx * (iy + mod(iz-1,Lz) * Ly);
        }

        u0_1p= Site[nn_ip].Psi[0].t - theta_1;
        u0_2p= Site[nn_ip].Psi[1].t - theta_2;

        u0_1m= -Site[nn_im].Psi[0].t + theta_1;
        u0_2m= -Site[nn_im].Psi[1].t + theta_2;

        for (n_1 = (-MCp.nMAX); n_1 < (MCp.nMAX+1); n_1++) {
            for (n_2 = (-MCp.nMAX); n_2 < (MCp.nMAX+1); n_2++) {

                u_1= u0_1p - C_TWO_PI*n_1;
                u_2= u0_2p - C_TWO_PI*n_2;

                local_S = 0.5 * my_beta * (Hp.rho * (u_1*u_1 + u_2*u_2) + Hp.nu*(u_1*u_2) );
                sum_n1n2_plus+= exp(-local_S);

                u_1= u0_1m - C_TWO_PI*n_1;
                u_2= u0_2m - C_TWO_PI*n_2;

                local_S = 0.5 * my_beta * (Hp.rho * (u_1*u_1 + u_2*u_2) + Hp.nu*(u_1*u_2) );
                sum_n1n2_minus+= exp(-local_S);
            }
        }

        boltz_h-=log(sum_n1n2_plus);
        boltz_h-=log(sum_n1n2_minus);

    }

    return boltz_h;

}




