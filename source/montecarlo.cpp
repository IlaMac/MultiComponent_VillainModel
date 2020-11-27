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
                for (alpha = 0; alpha < 3; alpha++) {
                    OldPsi = Site[i].Psi[alpha];
                    //t_localHtheta.tic();
                    oldE = HVillan_old(my_beta, Hp, MCp, Site);
                    //t_localHtheta.toc();
                    d_theta = rn::uniform_real_box(-MCp.lbox_theta, MCp.lbox_theta);
                    NewPsi.t = fmod(OldPsi.t + d_theta, C_TWO_PI);
                    NewPsi.r=1;
                    //t_localHtheta.tic();
                    newE = HVillan_new(newPsi, i, alpha, my_beta, Hp, MCp, Site)
                    //t_localHtheta.toc();
                    minus_deltaE =  (oldE - newE);
                    // h3 * (oldE - newE); when I will include properly the lattice spacing which in this case is taken to be just 1
                    if (minus_deltaE > 0) {
                        Site[i].Psi[alpha] = NewPsi;
                        acc_theta++;
                    } else {
                        rand = rn::uniform_real_box(0, 1);
                        //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
                        if (rand < exp(my_beta * minus_deltaE)) {
                            Site[i].Psi[alpha] = NewPsi;
                            acc_theta++;
                        }
                    }
                }
                if (Hp.e != 0) {
                    /**********VECTOR POTENTIAL UPDATE********/
                    for (vec = 0; vec < 3; vec++) {
                        //Update of A
                        OldA = Site[i].A[vec];
                        //t_localHA.tic();
                        oldE = local_HA(OldA, ix, iy, iz, vec, Hp, Site);
                        //t_localHA.toc();
                        d_A = rn::uniform_real_box(-MCp.lbox_A, MCp.lbox_A);
                        NewA = OldA + d_A;
                        //t_localHA.tic();
                        newE = local_HA(NewA, ix, iy, iz, vec, Hp, Site);
                        //t_localHA.toc();
                        minus_deltaE = h3 * (oldE - newE);
                        if (minus_deltaE > 0.) {
                            Site[i].A[vec] = NewA;
                            acc_A++;
                        } else {
                            rand = rn::uniform_real_box(0, 1);
                            //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
                            if (rand < exp(my_beta * minus_deltaE)) {
                                Site[i].A[vec] = NewA;
                                acc_A++;

                            }
                        }
                    }
                }
            }
        }
    }
//    t_localHA.print_measured_time();
 //   t_localHtheta.print_measured_time();
 //   t_localHpsi.print_measured_time();
    acc_theta=(double) acc_theta/(3*N);
    acc_A=(double) acc_A/(3*N);
    MCp.lbox_theta= MCp.lbox_theta*(0.5*(acc_theta/acc_rate)+0.5);
    MCp.lbox_A= MCp.lbox_A*(0.5*(acc_A/acc_rate)+0.5);
}

double HVillan_old(double my_beta, struct H_parameters &Hp, struct MC_parameters &MCp, struct Node* Site){

    unsigned int vec;
    unsigned int n_1, n_2;
    unsigned int i, ix, iy, iz, nn_i;
    double u_1, u_2;
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

                    u_1= Site[nn_i].Psi[1].t - Site[i].Psi[1].t;
                    u_2= Site[nn_i].Psi[2].t - Site[i].Psi[2].t;

                    for (n_1 = -MCp.INT_NMAX; n_1 < MCp.INT_NMAX; n_1++) {
                        for (n_2 = -MCp.INT_NMAX; n_2 < MCp.INT_NMAX; n_2++) {

                            u1-= C_TWO_PI*n_1;
                            u2-= C_TWO_PI*n_2;

                            local_S = 0.5 * my_beta * (Hp.rho * (u_1*u_1 + u_2*u_2) + Hp.nu*(u_1*u_2) );
                            sum_n1n2+= exp(-local_S);
                        }
                    }

                    boltz_h-=inv_beta*log(sum_n1n2);
                }
            }
        }
    }


    return boltz_h;

}

double HVillan_new(struct O2 Psi, unsigned int position, unsigned int alpha,  double my_beta, struct H_parameters &Hp, struct MC_parameters &MCp, struct Node* Site){

    unsigned int vec;
    unsigned int n_1, n_2;
    unsigned int i, ix, iy, iz, nn_i;
    double theta_1=0., theta_2=0.;
    double theta_nn_1=0., theta_nn_2=0.;
    double u_1, u_2;
    double local_S=0., sum_n1n2=0., boltz_h=0.;
    double inv_beta=1./my_beta;

    for(iz=0;iz<Lz;iz++){
        for(iy=0;iy<Ly;iy++){
            for(ix=0; ix<Lx; ix++){
                i = ix + Lx * (iy + iz * Ly);
                if(i == position){
                    if(alpha==0){
                        theta_1=Psi.t;
                        theta_1=Site[i].Psi[1].t;
                    }
                    else if(alpha==1){
                        theta_2=Psi.t;
                        theta_2=Site[i].Psi[2].t;
                    }
                }
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

                    if(nn_i == position){
                        if(alpha==0){
                            theta_nn_1=Psi.t;
                            theta_nn_1=Site[nn_i].Psi[1].t;
                        }
                        else if(alpha==1){
                            theta_nn_2=Psi.t;
                            theta_nn_2=Site[nn_i].Psi[2].t;
                        }
                    }

                    u_1= theta_nn_1 - theta_1;
                    u_2= theta_nn_2 - theta_2;

                    for (n_1 = -MCp.INT_NMAX; n_1 < MCp.INT_NMAX; n_1++) {
                        for (n_2 = -MCp.INT_NMAX; n_2 < MCp.INT_NMAX; n_2++) {

                            u1-= C_TWO_PI*n_1;
                            u2-= C_TWO_PI*n_2;

                            local_S = 0.5 * my_beta * (Hp.rho * (u_1*u_1 + u_2*u_2) + Hp.nu*(u_1*u_2) );
                            sum_n1n2+= exp(-local_S);
                        }
                    }

                    boltz_h-=inv_beta*log(sum_n1n2);
                }
            }
        }
    }

    return boltz_h;

}



double local_HA(double A, unsigned int ix, unsigned int iy, unsigned int iz,  unsigned int vec,  struct H_parameters &Hp, struct Node* Site){

    double h_Kinetic=0., h_B, h_AB=0., h_tot;
    double inv_h2=1./(Hp.h*Hp.h);
    double inv_h=1./(Hp.h);
    double J_alpha, J_beta, J_gamma;
    unsigned int alpha, i, l;
    double F2_A=0., F_A;
    double gauge_phase1;
    std::vector<struct O2> aux_field(3);
    unsigned int nn_ip, nn_im, nn_ipl, nn_iml;
    i=ix +Lx*(iy+Ly*iz);

    unsigned int ip=(ix == Lx-1 ? 0: ix+1);
    unsigned int ipx= ip+Lx*(iy+Ly*iz);
    unsigned int jp=(iy == Ly-1 ? 0: iy+1);
    unsigned int ipy= ix+(Lx*jp)+(Lx*Ly*iz);
    unsigned int kp=(iz == Lz-1 ? 0: iz+1);
    unsigned int ipz= ix+Lx*(iy+Ly*kp);

    unsigned int imx= (ix == 0 ? Lx-1: ix-1)+Lx*(iy+Ly*iz);
    unsigned int imy= ix+Lx*((iy == 0 ? Ly-1: iy-1)+Ly*iz);
    unsigned int imz= ix+Lx*(iy+Ly*(iz == 0 ? Lz-1: iz-1));

    if(vec==0){
        nn_ip=ipx;
        nn_im=imx;
    }else if(vec==1){
        nn_ip=ipy;
        nn_im=imy;
    }else if(vec==2) {
        nn_ip = ipz;
        nn_im = imz;
    }

    for(alpha=0; alpha<3; alpha++) {
        gauge_phase1=Site[nn_ip].Psi[alpha].t - Site[i].Psi[alpha].t + Hp.h*Hp.e*A;
        h_Kinetic -= inv_h2*cos(gauge_phase1);
    }

    //Compute the local Energy respect to a given component (alpha) of the vector potential A and a given spatial position (r=(ix, iy, iz))
    //We need to compute just the part of the Hamiltonian involving A

    //Kinetic= -(1/h²)*\sum_k=1,2,3 |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))
    for(alpha=0; alpha<3; alpha++) {
        gauge_phase1=Site[nn_ip].Psi[alpha].t - Site[i].Psi[alpha].t + Hp.h*Hp.e*A;
        h_Kinetic -= inv_h2*cos(gauge_phase1);
    }
    //Andreev-Bashkin term = \sum_beta!=alpha \sum_k=1,2,3 nu*(J^k_alpha - J^k_beta)^2;
    // with J^k_alpha= |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* sin(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r)))
    if(Hp.nu !=0 ) {

        J_alpha= inv_h * sin(Site[nn_ip].Psi[0].t - Site[i].Psi[0].t + Hp.h*Hp.e*A);
        J_beta= inv_h * sin(Site[nn_ip].Psi[1].t - Site[i].Psi[1].t + Hp.h*Hp.e*A);
        J_gamma= inv_h * sin(Site[nn_ip].Psi[2].t - Site[i].Psi[2].t + Hp.h*Hp.e*A);

        h_AB -= Hp.nu *( (J_alpha*J_beta)  + (J_alpha *J_gamma)+ (J_gamma*J_beta));
    }

    //All the plaquettes involving A_vec(i)
    for(l=0; l<3; l++){
        if(l==0){
            nn_ipl=ipx;
        }else if(l==1){
            nn_ipl=ipy;
        }else if(l==2) {
            nn_ipl = ipz;
        }
        if(l!= vec){
            F_A=(A + Site[nn_ip].A[l] - Site[nn_ipl].A[vec] - Site[i].A[l]);
            F2_A+=(F_A*F_A);
        }
    }
    for(l=0; l<3; l++){
        if(l==0){
            nn_iml=imx;
        }else if(l==1){
            nn_iml=imy;
        }else if(l==2) {
            nn_iml= imz;
        }
        if(l!= vec){
            F_A=(Site[nn_iml].A[vec]+ Site[nn(nn_iml, vec, 1)].A[l] - A - Site[nn_iml].A[l]);
            F2_A+=(F_A*F_A);
        }
    }

    h_B=(0.5*inv_h2)*F2_A;

    h_tot= h_Kinetic + h_AB+ h_B;
    return h_tot;
}

double F_2(double newA, unsigned int k, unsigned int ix, unsigned int iy, unsigned int iz, struct Node* Site){

    unsigned int l;
    unsigned int i;
    double F2_A=0., F_A;
    i=ix +Lx*(iy+Ly*iz);
    //All the plaquettes involving A_vec(i)
    for(l=0; l<3; l++){
        if(l!= k){
            F_A=(newA + Site[nn(i, k, 1)].A[l] - Site[nn(i, l, 1)].A[k] - Site[i].A[l]);
            F2_A+=(F_A*F_A);
        }
    }
    for(l=0; l<3; l++){
        if(l!= k){
            F_A=(Site[nn(i, l, -1)].A[k]+ Site[nn(nn(i, l, -1), k, 1)].A[l] - newA - Site[nn(i, l, -1)].A[l]);
            F2_A+=(F_A*F_A);
        }
    }

    return F2_A;
}

