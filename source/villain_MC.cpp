//
// Created by ilaria on 2020-12-17.
//

#include "villain_MC.h"


void growCluster(int i, int* clusterSpin, const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta, struct Villain &vil){

    int new_Psi[NC]={0};
    int old_Psi[NC]={0};
    int ix, iy, iz;
    int ip, im, vec, alpha;
    int beta;
    int arg_F_new[NC]={0};
    int arg_B_new[NC]={0}; /*Forward and Backward updated phases*/
    int arg_F_old[NC]={0};
    int arg_B_old[NC]={0}; /*Forward and Backward updated phases*/
    //addProbability for the link i,j has the form: 1- exp(-beta \Delta E_ij) with \Delta E_ij
    double rand=0., addProbability;
    double dE=0;

    ix= i%Lx;
    iy= (i/Lx)%Ly;
    iz= i/(Lx*Ly);


    //phase in i added to the cluster and flipped
    clusterSpin[i]=1;
    int phase_diff = Site[i].Psi[0] - Site[i].Psi[1];

    old_Psi[0]= Site[i].Psi[0];
    old_Psi[1]= Site[i].Psi[1];

    new_Psi[0] = Site[i].Psi[1];
    new_Psi[1] = Site[i].Psi[0];

    Site[i].Psi[0]= new_Psi[0];
    Site[i].Psi[1]= new_Psi[1];

    //Check of the neighbours of i
    // if the neighbor spin does not belong to the cluster, but it has the preconditions to be added, then try to add it to the cluster
    for (vec = 0; vec < 3; vec++) {
        if (vec == 0) {
            ip = mod(ix + 1, Lx) + Lx * (iy + iz * Ly);
            im = mod(ix - 1, Lx) + Lx * (iy + iz * Ly);
        }
        if (vec == 1) {
            ip = ix + Lx * (mod(iy + 1, Ly) + iz * Ly);
            im = ix + Lx * (mod(iy - 1, Ly) + iz * Ly);
        }
        if (vec == 2) {
            ip = ix + Lx * (iy + mod(iz + 1, Lz) * Ly);
            im = ix + Lx * (iy + mod(iz - 1, Lz) * Ly);
        }
        if (clusterSpin[im]==0){
            for(alpha=0; alpha<NC; alpha++){
                arg_B_new[alpha]  = arg( (new_Psi[alpha] - Site[im].Psi[alpha]) - inv_dp*Hp.e*Site[im].A[vec], MaxP);
                arg_B_old[alpha] = arg( (old_Psi[alpha]  - Site[im].Psi[alpha]) - inv_dp*Hp.e*Site[im].A[vec], MaxP);
            }
            //with this sweep I am not touching the Josephson term.
            dE = (  vil.potential.at(OFFSET_POT + arg_B_new[0] + MaxP * arg_B_new[1]) -vil.potential.at(OFFSET_POT + arg_B_old[0] + MaxP * arg_B_old[1]));
            addProbability=1-exp(-my_beta*dE);
            rand= rn::uniform_real_box(0, 1);;
            if (rand < addProbability){
                growCluster(im, clusterSpin, Site, MCp, Hp, my_beta, vil);}
        }
        if (clusterSpin[ip]==0){
            for(alpha=0; alpha<NC; alpha++){
                arg_F_new[alpha] = arg( (Site[ip].Psi[alpha] - new_Psi[alpha]) - inv_dp*Hp.e*Site[i].A[vec], MaxP);
                arg_F_old[alpha] = arg( (Site[ip].Psi[alpha] - old_Psi[alpha]) - inv_dp*Hp.e*Site[i].A[vec], MaxP);
            }
            //with this sweep I am not touching the Josephson term.
            dE = (  vil.potential.at(OFFSET_POT + arg_F_new[0] + MaxP * arg_F_new[1]) - vil.potential.at(OFFSET_POT + arg_F_old[0] + MaxP * arg_F_old[1]));
            addProbability=1-exp(-my_beta*dE);

            rand= rn::uniform_real_box(0, 1);;
            if (rand < addProbability){
                growCluster(ip, clusterSpin, Site, MCp, Hp, my_beta, vil);}
        }
    }

    return;
}

void wolff(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta, struct Villain &vil){

    double rand=0;
    int i, iseed, count=0;
    int clusterSpin[N]={0};
    double  phase_diff, phase_diff_seed;

    /*choose randomly a site of the lattice*/
    iseed = rn::uniform_integer_box(0, N-1);
    phase_diff_seed=dp*(Site.at(iseed).Psi[1] - Site.at(iseed).Psi[0]);
    while(phase_diff_seed > M_PI){
        phase_diff_seed-= 2*M_PI;}
    while(phase_diff_seed<=-M_PI){
        phase_diff_seed+=2*M_PI;}

//    for(int iz=0; iz<Lz; iz++){
//        for(int iy=0; iy<Ly; iy++){
//            for(int ix=0; ix<Lx; ix++){
//                i=ix+ Lx*(iy +iz*Ly);
//                phase_diff=dp*(Site.at(i).Psi[1] - Site.at(i).Psi[0]);
//                while(phase_diff > M_PI){
//                    phase_diff-= 2*M_PI;}
//                while(phase_diff<=-M_PI){
//                    phase_diff+=2*M_PI;}
//                if( (phase_diff_seed)*(phase_diff)<0){
//                    clusterSpin[i]=-1;
//                }
//            }
//        }
//    }

    growCluster(iseed, clusterSpin, Site, MCp, Hp, my_beta, vil);

    for(int iz=0; iz<Lz; iz++){
        for(int iy=0; iy<Ly; iy++){
            for(int ix=0; ix<Lx; ix++){
                i=ix+ Lx*(iy +iz*Ly);
                if( clusterSpin[i]==1){count++;}
            }
        }
    }

//    std::cout<< "WOLFF size: "<< count<<std::endl;

    return;
}



void metropolis_villain(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta, struct Villain &vil){

    int ix, iy, iz;
    int ipx, ipy, ipz;
    int ip, im, alpha, vec, i;
    int n_var;
    int new_int_phase[NC]={0};
    int arg_F_new[NC][3]={{0}};
    int arg_B_new[NC][3]={{0}}; /*Forward and Backward updated phases*/
    int arg_F_old[NC][3]={{0}};
    int arg_B_old[NC][3]={{0}}; /*Forward and Backward updated phases*/
    double acc_rate=0.5, acc_theta=0., acc_locked=0., rand;
    double dE;
    double dE_A;
    double OldA, NewA, dA, acc_A=0.;
    double F_A_new, B_A_new, F_A_old, B_A_old;
    int l, nn_ipl, nn_iml, nn_ipvec_iml;
    double curl2_A_new, curl2_A_old;

    class_tic_toc t_localHtheta(true,5,"local_Htheta");

    for (iz= 0; iz < Lz; iz++) {
        for (iy = 0; iy < Ly; iy++) {
            for (ix = 0; ix < Lx; ix++) {
                i = ix + Lx * (iy + iz * Ly);

//                std::cout<< "site: "<< i <<"cos phase_diff: "<< cos(2*dp*(Site[i].Psi[0] - Site[i].Psi[1])) << "phase diff/pi: "<<dp*(Site[i].Psi[0] - Site[i].Psi[1])/M_PI << std::endl;
                /*******INDIVIDUAL PHASES ONLY UPDATE**************/
                for (alpha = 0; alpha < NC; alpha++) {
                    n_var = rn::uniform_integer_box(1, MCp.lbox);
                    //Try to extract a new value of the phase directly in the interval [-(MaxP-1)/2; (MaxP-1/2]
                    rand = rn::uniform_real_box(0, 1);
                    if(rand<0.5){
                        new_int_phase[alpha]= arg(Site[i].Psi[alpha] + n_var, MaxP);}
                    else{
                        new_int_phase[alpha]= arg(Site[i].Psi[alpha] - n_var, MaxP);
                    }
                    //std::cout<< "Component: "<< alpha << " Old phase: "<<  Site[i].Psi[alpha] << " New phase: " << new_int_phase[alpha] << " increment: " << n_var << std::endl;
                    for (vec = 0; vec < 3; vec++) {
                        if (vec == 0) {
                            ip = mod(ix + 1, Lx) + Lx * (iy + iz * Ly);
                            im = mod(ix - 1, Lx) + Lx * (iy + iz * Ly);
                        }
                        if (vec == 1) {
                            ip = ix + Lx * (mod(iy + 1, Ly) + iz * Ly);
                            im = ix + Lx * (mod(iy - 1, Ly) + iz * Ly);
                        }
                        if (vec == 2) {
                            ip = ix + Lx * (iy + mod(iz + 1, Lz) * Ly);
                            im = ix + Lx * (iy + mod(iz - 1, Lz) * Ly);
                        }
                        //std::cout<< "ix: "<< ix << " iy: "<< iy <<" iz: "<< iz << " vec= "<< vec << " ip: " << ip << " im: "<< im << std::endl;
                         arg_F_new[alpha][vec] = arg( (Site[ip].Psi[alpha] - new_int_phase[alpha]) - inv_dp*Hp.e*Site[i].A[vec], MaxP);
                         arg_B_new[alpha][vec] = arg( (new_int_phase[alpha] - Site[im].Psi[alpha]) - inv_dp*Hp.e*Site[im].A[vec], MaxP);
                         arg_F_old[alpha][vec] = arg( (Site[ip].Psi[alpha] - Site[i].Psi[alpha]) - inv_dp*Hp.e*Site[i].A[vec], MaxP);
                         arg_B_old[alpha][vec] = arg( (Site[i].Psi[alpha] - Site[im].Psi[alpha]) - inv_dp*Hp.e*Site[im].A[vec], MaxP);
                    }
                }
                /******************Specific for the two component case******************/

                /*******************Phase updating phase componet 1*******************/
                dE=0;
                for (vec = 0; vec < DIM; vec++) {
                   dE += (vil.potential.at(OFFSET_POT + arg_B_new[0][vec] + MaxP * arg_B_old[1][vec])
                         +vil.potential.at(OFFSET_POT + arg_F_new[0][vec] + MaxP * arg_F_old[1][vec])
                         -vil.potential.at(OFFSET_POT + arg_B_old[0][vec] + MaxP * arg_B_old[1][vec])
                         -vil.potential.at(OFFSET_POT + arg_F_old[0][vec] + MaxP * arg_F_old[1][vec]));
                }
                rand = rn::uniform_real_box(0, 1);
                dE+=( (Hp.eta1*cos(dp*(new_int_phase[0] - Site[i].Psi[1])) + Hp.eta2*cos(2*dp*(new_int_phase[0] - Site[i].Psi[1])) )
                        -(Hp.eta1*cos(dp*(Site[i].Psi[0] - Site[i].Psi[1])) + Hp.eta2*cos(2*dp*(Site[i].Psi[0] - Site[i].Psi[1]))) );

                //Boltzmann weight: exp(-\beta dH) dH= 1/beta dE
                //std::cout<< rand << " exp mybetadE: "<< exp(-my_beta*dE)<< std::endl;
                //std::cout<< -log(rand) << " mybetadE: "<< my_beta*dE<< std::endl;
                if(dE<0){
                    Site[i].Psi[0] = new_int_phase[0];
                    acc_theta++;
                    for(vec=0; vec<3;vec++){
                        arg_B_old[0][vec]=arg_B_new[0][vec];
                        arg_F_old[0][vec]=arg_F_new[0][vec];
                    }
                }else{
                    if (rand <exp(-my_beta*dE)) {
                        Site[i].Psi[0] = new_int_phase[0];
                        acc_theta++;
                        for(vec=0; vec<3;vec++){
                            arg_B_old[0][vec]=arg_B_new[0][vec];
                            arg_F_old[0][vec]=arg_F_new[0][vec];
                        }
                    }
                }

                /*******************Phase updating phase componet 2*******************/
                dE=0;
                for (vec = 0; vec < DIM; vec++) {
                   dE += (vil.potential.at(OFFSET_POT + arg_B_old[0][vec] + MaxP * arg_B_new[1][vec])
                         +vil.potential.at(OFFSET_POT + arg_F_old[0][vec] + MaxP * arg_F_new[1][vec])
                         -vil.potential.at(OFFSET_POT + arg_B_old[0][vec] + MaxP * arg_B_old[1][vec])
                         -vil.potential.at(OFFSET_POT + arg_F_old[0][vec] + MaxP * arg_F_old[1][vec]));
                }
                dE+=( (Hp.eta1*cos(dp*(new_int_phase[1] - Site[i].Psi[0])) + Hp.eta2*cos(2*dp*(new_int_phase[1] - Site[i].Psi[0])))
                        -(Hp.eta1*cos(dp*(Site[i].Psi[0] - Site[i].Psi[1])) + Hp.eta2*cos(2*dp*(Site[i].Psi[0] - Site[i].Psi[1]))) );

                if(dE<0){
                    Site[i].Psi[1] = new_int_phase[1];
                    acc_theta++;
                }else{
                    rand = rn::uniform_real_box(0, 1);
                    //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
                    if (rand < exp(-my_beta*dE)) {
                        Site[i].Psi[1] = new_int_phase[1];
                        acc_theta++;
                    }
                }


                /*******************PHASE UPDATING OF BOTH PHASES BY THE SAME PHASE SHIFT*******************/
                /*This kind of updates does not change the local Josephson energy since the phase difference is not modified*/

                n_var = rn::uniform_integer_box(1, MCp.lbox);
                rand = rn::uniform_real_box(0, 1);
                for (alpha = 0; alpha < NC; alpha++) {
                    if(rand<0.5){
                        new_int_phase[alpha]= arg(Site[i].Psi[alpha] + n_var, MaxP);}
                    else{
                        new_int_phase[alpha]= arg(Site[i].Psi[alpha] - n_var, MaxP);
                    }

                    for (vec = 0; vec < 3; vec++) {
                        if (vec == 0) {
                            ip = mod(ix + 1, Lx) + Lx * (iy + iz * Ly);
                            im = mod(ix - 1, Lx) + Lx * (iy + iz * Ly);
                        }
                        if (vec == 1) {
                            ip = ix + Lx * (mod(iy + 1, Ly) + iz * Ly);
                            im = ix + Lx * (mod(iy - 1, Ly) + iz * Ly);
                        }
                        if (vec == 2) {
                            ip = ix + Lx * (iy + mod(iz + 1, Lz) * Ly);
                            im = ix + Lx * (iy + mod(iz - 1, Lz) * Ly);
                        }
                        //std::cout<< "ix: "<< ix << " iy: "<< iy <<" iz: "<< iz << " vec= "<< vec << " ip: " << ip << " im: "<< im << std::endl;
                        arg_F_new[alpha][vec] = arg( (Site[ip].Psi[alpha] - new_int_phase[alpha]) - inv_dp*Hp.e*Site[i].A[vec], MaxP);
                        arg_B_new[alpha][vec] = arg( (new_int_phase[alpha] - Site[im].Psi[alpha]) - inv_dp*Hp.e*Site[im].A[vec], MaxP);
                        arg_F_old[alpha][vec] = arg( (Site[ip].Psi[alpha] - Site[i].Psi[alpha]) - inv_dp*Hp.e*Site[i].A[vec], MaxP);
                        arg_B_old[alpha][vec] = arg( (Site[i].Psi[alpha] - Site[im].Psi[alpha]) - inv_dp*Hp.e*Site[im].A[vec], MaxP);
                    }
                }
                dE=0;
                for (vec = 0; vec < DIM; vec++) {
                    dE += (vil.potential.at(OFFSET_POT + arg_B_new[0][vec] + MaxP * arg_B_new[1][vec])
                            +vil.potential.at(OFFSET_POT + arg_F_new[0][vec] + MaxP * arg_F_new[1][vec])
                            -vil.potential.at(OFFSET_POT + arg_B_old[0][vec] + MaxP * arg_B_old[1][vec])
                            -vil.potential.at(OFFSET_POT + arg_F_old[0][vec] + MaxP * arg_F_old[1][vec]));
                }
//                std::cout<<"dE prima: "<< dE<<std::endl;
//                dE+=( (Hp.eta1*cos(dp*(new_int_phase[0] - new_int_phase[1])) + Hp.eta2*cos(2*dp*(new_int_phase[0] - new_int_phase[1])) )
//                        -(Hp.eta1*cos(dp*(Site[i].Psi[0] - Site[i].Psi[1])) + Hp.eta2*cos(2*dp*(Site[i].Psi[0] - Site[i].Psi[1]))) );
//                std::cout<<"dE dopo: "<< dE<<std::endl;

                if(dE<0){
                    Site[i].Psi[0] = new_int_phase[0];
                    Site[i].Psi[1] = new_int_phase[1];
                    acc_locked++;
                }else{
                    //Boltzmann weight: exp(-\beta dH) dH= 1/beta dE
                    rand = rn::uniform_real_box(0, 1);
                    if (rand <exp(-my_beta*dE)) {
                        Site[i].Psi[0] = new_int_phase[0];
                        Site[i].Psi[1] = new_int_phase[1];
                        acc_locked++;
                    }
                }

                /*******************VECTOR POTENTIAL UPDATE*******************/
                if(Hp.e!=0){
                    for (vec=0; vec<DIM; vec++){
                        OldA= Site[i].A[vec];
                        dA = rn::uniform_real_box(-MCp.lbox_A, MCp.lbox_A);
                        NewA = OldA + dA;
                        ipx=ix;
                        ipy=iy;
                        ipz=iz;
                        if (vec == 0) {
                            ipx=mod(ix + 1, Lx);
                        }
                        if (vec == 1) {
                            ipy=mod(iy + 1, Ly);
                        }
                        if (vec == 2) {
                            ipz = mod(iz + 1, Lz);
                        }
                        ip = ipx + Lx * (ipy + ipz * Ly);
                        for(alpha=0; alpha<NC; alpha++){
                            arg_F_new[alpha][vec] = arg( (Site[ip].Psi[alpha] - Site[i].Psi[alpha]) - inv_dp*Hp.e*NewA, MaxP);
                            arg_F_old[alpha][vec] = arg( (Site[ip].Psi[alpha] - Site[i].Psi[alpha]) - inv_dp*Hp.e*OldA, MaxP);
                        }
                        dE_A= vil.potential.at(OFFSET_POT + arg_F_new[0][vec] + MaxP * arg_F_new[1][vec]) - vil.potential.at(OFFSET_POT + arg_F_old[0][vec] + MaxP * arg_F_old[1][vec]);
                        //curl A
                        curl2_A_new=0.;
                        curl2_A_old=0.;
                        //All the plaquettes involving A_vec(i)
                        for(l=0; l<DIM; l++){
                            if(l!= vec){
                                if(l==0){
                                    nn_ipl= mod(ix + 1, Lx) + Lx * (iy + iz * Ly);
                                    nn_iml= mod(ix - 1, Lx) + Lx * (iy + iz * Ly);
                                    nn_ipvec_iml= mod(ipx - 1, Lx) + Lx * (ipy + ipz * Ly);
                                }else if(l==1){
                                    nn_ipl= ix + Lx * (mod(iy + 1, Ly) + iz * Ly);
                                    nn_iml= ix + Lx * (mod(iy - 1, Ly) + iz * Ly);
                                    nn_ipvec_iml= ipx + Lx * (mod(ipy - 1, Ly) + ipz * Ly);
                                }else if(l==2) {
                                    nn_ipl= ix + Lx * (iy + mod(iz + 1, Lz) * Ly);
                                    nn_iml= ix + Lx * (iy + mod(iz - 1, Lz) * Ly);
                                    nn_ipvec_iml= ipx + Lx * (ipy + mod(ipz - 1, Lz) * Ly);
                                }
                                //plaquette to the right (F="forward") side of newA
                                F_A_new=(NewA + Site[ip].A[l] - Site[nn_ipl].A[vec] - Site[i].A[l]);
                                curl2_A_new+=0.5*(F_A_new*F_A_new);
                                //plaquette to the left (B="backrward") side of newA
                                B_A_new=(Site[nn_iml].A[vec]+ Site[nn_ipvec_iml].A[l] - NewA - Site[nn_iml].A[l]);
                                curl2_A_new+=0.5*(B_A_new*B_A_new);

                                F_A_old=(OldA + Site[ip].A[l] - Site[nn_ipl].A[vec] - Site[i].A[l]);
                                curl2_A_old+=0.5*(F_A_old*F_A_old);
                                B_A_old=(Site[nn_iml].A[vec]+ Site[nn_ipvec_iml].A[l] - OldA - Site[nn_iml].A[l]);
                                curl2_A_old+=0.5*(B_A_old*B_A_old);
                            }
                        }
                        dE_A+= curl2_A_new -curl2_A_old;
                        if(dE_A<0){
                            Site[i].A[vec] = NewA;
                            acc_A++;
                        }else{
                            rand = rn::uniform_real_box(0, 1);
                            //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
                            if (rand < exp(-my_beta*dE_A)) {
                                Site[i].A[vec] = NewA;
                                acc_A++;
                            }
                        }
                    }
                }
            }
        }
    }

    acc_theta=(double) acc_theta/(N*NC);
    acc_locked=(double) acc_locked/(N);
    acc_A=(double) acc_A/(DIM*N);

    if(acc_theta/acc_rate>1){MCp.lbox+=1;}
    if( (acc_theta/acc_rate<1) and (MCp.lbox>1)){ MCp.lbox-=1;}

    if(acc_locked/acc_rate>1){MCp.lbox_coupled+=1;}
    if( (acc_locked/acc_rate<1) and (MCp.lbox_coupled>1)){ MCp.lbox_coupled-=1;}
    MCp.lbox_A= MCp.lbox_A*((0.5*acc_A/acc_rate)+0.5);

//    MPI_Barrier(MPI_COMM_WORLD);
//    std::cout<<" beta: "<< my_beta << std::endl;
//    std::cout<<" acc_theta: "<< acc_theta << " MCp_lbox: " << MCp.lbox << std::endl;
//    std::cout<<" acc_locked: "<< acc_locked << " MCp_lbox_coupled: " << MCp.lbox_coupled << std::endl;
//    std::cout<<" acc_theta1: "<< acc_theta1 << " acc_theta2: "<< acc_theta2 <<" acc_locked: " <<  acc_locked << std::endl;

//    //New Metropolis sweeps consisting in the reflection of one phase respect to the other so to switch between \phi_12 = \pi/2 and \phi_12 = -\pi/2
//    for (iz= 0; iz < Lz; iz++) {
//        for (iy = 0; iy < Ly; iy++) {
//            for (ix = 0; ix < Lx; ix++) {
//                i = ix + Lx * (iy + iz * Ly);
//                new_int_phase[0]= Site[i].Psi[0] ;
//                new_int_phase[1]= Site[i].Psi[1] ;
//                alpha=rn::uniform_integer_box(0,1);
//                int phase_diff= Site[i].Psi[0] - Site[i].Psi[1];
//                // if alpha=0: phi_1 -> phi_1 -2*phase_diff; if alpha=1: phi_2 -> phi_2 +2*phase_diff;
//                new_int_phase[alpha]= arg(Site[i].Psi[alpha] + 2*(2*alpha -1)*phase_diff, MaxP);
//                int new_phase_diff=new_int_phase[0] -new_int_phase[1];
//
//                for (vec = 0; vec < 3; vec++) {
//                    if (vec == 0) {
//                        ip = mod(ix + 1, Lx) + Lx * (iy + iz * Ly);
//                        im = mod(ix - 1, Lx) + Lx * (iy + iz * Ly);
//                    }
//                    if (vec == 1) {
//                        ip = ix + Lx * (mod(iy + 1, Ly) + iz * Ly);
//                        im = ix + Lx * (mod(iy - 1, Ly) + iz * Ly);
//                    }
//                    if (vec == 2) {
//                        ip = ix + Lx * (iy + mod(iz + 1, Lz) * Ly);
//                        im = ix + Lx * (iy + mod(iz - 1, Lz) * Ly);
//                    }
//                    //std::cout<< "ix: "<< ix << " iy: "<< iy <<" iz: "<< iz << " vec= "<< vec << " ip: " << ip << " im: "<< im << std::endl;
//
//                    arg_F_new[0][vec] = arg( (Site[ip].Psi[0] - new_int_phase[0]) - inv_dp*Hp.e*Site[i].A[vec], MaxP);
//                    arg_B_new[0][vec] = arg( (new_int_phase[0] - Site[im].Psi[0]) - inv_dp*Hp.e*Site[im].A[vec], MaxP);
//                    arg_F_old[0][vec] = arg( (Site[ip].Psi[0] - Site[i].Psi[0]) - inv_dp*Hp.e*Site[i].A[vec], MaxP);
//                    arg_B_old[0][vec] = arg( (Site[i].Psi[0] - Site[im].Psi[0]) - inv_dp*Hp.e*Site[im].A[vec], MaxP);
//
//                    arg_F_new[1][vec] = arg( (Site[ip].Psi[1] - new_int_phase[1]) - inv_dp*Hp.e*Site[i].A[vec], MaxP);
//                    arg_B_new[1][vec] = arg( (new_int_phase[1] - Site[im].Psi[1]) - inv_dp*Hp.e*Site[im].A[vec], MaxP);
//                    arg_F_old[1][vec] = arg( (Site[ip].Psi[1] - Site[i].Psi[1]) - inv_dp*Hp.e*Site[i].A[vec], MaxP);
//                    arg_B_old[1][vec] = arg( (Site[i].Psi[1] - Site[im].Psi[1]) - inv_dp*Hp.e*Site[im].A[vec], MaxP);
//                }
//                dE= 0;
//                //with this sweep I am not touching the Josephson term.
//                for (vec = 0; vec < DIM; vec++) {
//                    dE += (vil.potential.at(OFFSET_POT + arg_B_new[0][vec] + MaxP * arg_B_new[1][vec])
//                            +vil.potential.at(OFFSET_POT + arg_F_new[0][vec] + MaxP * arg_F_new[1][vec])
//                            -vil.potential.at(OFFSET_POT + arg_B_old[0][vec] + MaxP * arg_B_old[1][vec])
//                            -vil.potential.at(OFFSET_POT + arg_F_old[0][vec] + MaxP * arg_F_old[1][vec]));
//                }
//                if(dE<0){
//                    Site[i].Psi[0] = new_int_phase[0];
//                    Site[i].Psi[1] = new_int_phase[1];
//                }else{
//                    //Boltzmann weight: exp(-\beta dH) dH= 1/beta dE
//                    rand = rn::uniform_real_box(0, 1);
//                    if (rand <exp(-my_beta*dE)) {
//                        Site[i].Psi[0] = new_int_phase[0];
//                        Site[i].Psi[1] = new_int_phase[1];
//                    }
//                }
//        }
//    }
//}

}



int arg(int x, int Max){

    while(x < -(Max -1)/2){
        x+=(Max-1);
    }
    while(x > (Max -1)/2){
        x-=(Max-1);
    }
    return x;
}