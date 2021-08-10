//
// Created by ilaria on 2019-11-16.
//

#include "measures.h"

void gauge_potential(struct Measures &mis, const std::vector<Node> &Site, struct H_parameters &Hp){
    int vec, l;
    int nn_i, nn_il;

    for(int iz=0;iz<Lz;iz++){
        for(int iy=0;iy<Ly;iy++){
            for(int ix=0; ix<Lx; ix++) {
                int i = ix + Lx * (iy + iz * Ly);
                for(vec=0; vec<3; vec++) {
                    if(vec==0){
                        nn_i= mod(ix+1, Lx) + Lx * (iy + iz * Ly);
                    }
                    if(vec==1){
                        nn_i= ix + Lx * (mod(iy+1,Ly) + iz * Ly);
                    }
                    if(vec==2){
                        nn_i= ix + Lx * (iy + mod(iz+1,Lz) * Ly);
                    }
                    for(l=vec; l<DIM; l++){
                        if(l==0){
                            nn_il= mod(ix+1, Lx) + Lx * (iy + iz * Ly);
                        }
                        if(l==1){
                            nn_il= ix + Lx * (mod(iy+1,Ly) + iz * Ly);
                        }
                        if(l==2){
                            nn_il= ix + Lx * (iy + mod(iz+1,Lz) * Ly);
                        }
                        //F_{vec, l}= A_vec(r_i) + A_l(ri+vec) - A_vec(r_i+l) - A_l(ri)
                        double F_A = (Site[i].A[vec] + Site[nn_i].A[l] - Site[nn_il].A[vec] -Site[i].A[l]);
                        mis.E_gp += (0.5 * (F_A * F_A));
                    }
                }
            }
        }
    }
}

void josephson_potential(struct Measures &mis, const std::vector<Node> &Site, struct H_parameters &Hp){
    for(int iz=0;iz<Lz;iz++){
        for(int iy=0;iy<Ly;iy++){
            for(int ix=0; ix<Lx; ix++) {
                int i = ix + Lx * (iy + iz * Ly);
                mis.E_jp+=Hp.eta1*cos(dp*(Site[i].Psi[0] - Site[i].Psi[1])) + Hp.eta2*cos(2*dp*(Site[i].Psi[0] - Site[i].Psi[1]));
            }
        }
    }

}

void energies(struct Measures &mis, struct Villain &vil, double &E_betanp, double &E_betanm, const std::vector<Node> &Site, struct H_parameters &Hp){

    double E_betanm_temp=0.;
    double E_betanp_temp=0.;

    for(int iz=0;iz<Lz;iz++){
        for(int iy=0;iy<Ly;iy++){
            for(int ix=0; ix<Lx; ix++) {
                int i = ix + Lx * (iy + iz * Ly);
                for(int vec=0; vec<3; vec++) {
                    int ipx=ix;
                    int ipy=iy;
                    int ipz=iz;
                    if (vec == 0) {
                        ipx=mod(ix + 1, Lx);
                    }
                    if (vec == 1) {
                        ipy=mod(iy + 1, Ly);
                    }
                    if (vec == 2) {
                        ipz = mod(iz + 1, Lz);
                    }
                   int nn_i = ipx + Lx * (ipy + ipz * Ly);

                    int arg_1 = arg(Site[nn_i].Psi[0] - Site.at(i).Psi[0] - inv_dp*Hp.e * Site[i].A[vec], MaxP);
                    int arg_2 = arg(Site[nn_i].Psi[1] - Site.at(i).Psi[1] - inv_dp*Hp.e * Site[i].A[vec], MaxP);
                    mis.U+= vil.upotential.at(OFFSET_POT + arg_1 +MaxP*arg_2);
                    mis.E+= vil.potential.at(OFFSET_POT + arg_1 +MaxP*arg_2) ;
                    E_betanm_temp+= vil.potential_bminus.at(OFFSET_POT + arg_1 +MaxP*arg_2);
                    E_betanp_temp+= vil.potential_bplus.at(OFFSET_POT + arg_1 +MaxP*arg_2);
                }
            }
        }
    }


    /****Gauge term***/
    if(Hp.e!=0) {
        gauge_potential(mis, Site, Hp);
    }
    /****Josephson term***/
    if((Hp.eta1!=0) or (Hp.eta2!=0)) {
        josephson_potential(mis, Site, Hp);
    }

    mis.U+=mis.E_gp + mis.E_jp;
    mis.E+=mis.E_gp + mis.E_jp;

    E_betanp=E_betanp_temp + mis.E_jp + mis.E_gp;
    E_betanm=E_betanm_temp + mis.E_jp + mis.E_gp;

}

void helicity_modulus(struct Measures &mis, struct Villain &vil, const std::vector<Node> &Site, struct H_parameters &Hp){

    //I compute the helicity modulus only along one direction x
    double d1=0., d2=0., d11=0., d22=0., d12=0.;

    for(int iz=0;iz<Lz;iz++){
        for(int iy=0;iy<Ly;iy++){
            for(int ix=0; ix<Lx; ix++) {
                int i = ix + Lx * (iy + iz * Ly);
                int nn_i = mod(ix + 1, Lx) + Lx * (iy + iz * Ly);

                int arg_1 = arg(Site[nn_i].Psi[0] - Site.at(i).Psi[0] - inv_dp* Hp.e * Site[i].A[0], MaxP);
                int arg_2 = arg(Site[nn_i].Psi[1] - Site.at(i).Psi[1] - inv_dp* Hp.e * Site[i].A[0], MaxP);

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

////DESIGNED FOR 2 COMPONENTs
void magnetization(struct Measures &mis, const std::vector<Node> &Site){
    //The Ising parameter m(x,y)=+/-1 indicates the chirality of the two phases. If the phase difference is + \pi/2 then m=1; otherwise  m=-1.
    int ix, iy, iz, i;
    double phi_shifted_1;

    for(iz=0; iz<Lz;iz++) {
        for (iy = 0; iy < Ly; iy++) {
            for (ix = 0; ix < Lx; ix++) {
                i=ix +Lx*(iy+Ly*iz);

                phi_shifted_1= dp*(Site.at(i).Psi[1] - Site.at(i).Psi[0]);
                while(phi_shifted_1 > M_PI){
                        phi_shifted_1-= 2*M_PI;}
                while(phi_shifted_1<=-M_PI){
                        phi_shifted_1+=2*M_PI;}

                if(phi_shifted_1>0){
                    mis.m+=1;
                }else if(phi_shifted_1<0){
                    mis.m+=(-1);
                }
            }
        }
    }
    mis.m=mis.m/N;
}

void magnetization_singlephase(struct Measures &mis, const std::vector<Node> &Site){
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

void dual_stiffness(struct Measures &mis, struct H_parameters &Hp, const std::vector<Node> &Site){
    double qx_min=C_TWO_PI/(Lx);
    double invNorm= 1./((C_TWO_PI)*(C_TWO_PI)*N);
    double Re_rhoz=0.;
    double Im_rhoz=0.;
    double inv_h=1./(Hp.h);

    for(int ix=0; ix<Lx;ix++){
        for(int iy=0; iy<Ly;iy++){
            for(int iz=0; iz<Lz;iz++){
                int i = ix + Lx * (iy + iz * Ly);
                int ipx = mod(ix + 1, Lx) + Lx * (iy + iz * Ly);
                int ipy = ix + Lx * (mod(iy + 1, Ly) + iz * Ly);
                double Dx_Ay=(Site[ipx].A[1]- Site[i].A[1])*inv_h;
                double Dy_Ax=(Site[ipy].A[0]- Site[i].A[0])*inv_h;
                Re_rhoz+=(cos((double)qx_min*ix)*(Dx_Ay -Dy_Ax));
                Im_rhoz+=(sin((double)qx_min*ix)*(Dx_Ay -Dy_Ax));
            }
        }
    }
    mis.d_rhoz=invNorm*((Re_rhoz*Re_rhoz) +(Im_rhoz*Im_rhoz));
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
            fwrite(s.A.data(), sizeof(double), DIM, fA);
	}
        fclose(fA);
    }

}


//
//void u_internal_energy(struct Measures &mis, struct Villain &vil, const std::vector<Node> &Site, struct H_parameters &Hp) {
//
//    int vec;
//    int arg_1, arg_2;
//    int i, ix, iy, iz, ipx, ipy, ipz, nn_i;
//
//    for(iz=0;iz<Lz;iz++){
//        for(iy=0;iy<Ly;iy++){
//            for(ix=0; ix<Lx; ix++) {
//                i = ix + Lx * (iy + iz * Ly);
//                for(vec=0; vec<3; vec++) {
//                    ipx=ix;
//                    ipy=iy;
//                    ipz=iz;
//                    if (vec == 0) {
//                        ipx=mod(ix + 1, Lx);
//                    }
//                    if (vec == 1) {
//                        ipy=mod(iy + 1, Ly);
//                    }
//                    if (vec == 2) {
//                        ipz = mod(iz + 1, Lz);
//                    }
//                    nn_i = ipx + Lx * (ipy + ipz * Ly);
//
//                    arg_1 = arg(Site[nn_i].Psi[0] - Site.at(i).Psi[0] - inv_dp*Hp.e * Site[i].A[vec], MaxP);
//                    arg_2 = arg(Site[nn_i].Psi[1] - Site.at(i).Psi[1] - inv_dp*Hp.e * Site[i].A[vec], MaxP);
//                    mis.U+= vil.upotential.at(OFFSET_POT + arg_1 +MaxP*arg_2);
//
//                }
//            }
//        }
//    }
//
//    /****Gauge term***/
//    mis.U+=mis.E_gp;
//    /****Josephson term***/
//    mis.U+=mis.E_jp;
//}
//
//void energy_nn(struct Measures &mis, struct Villain &vil, double &E_betanp, double &E_betanm, const std::vector<Node> &Site, struct H_parameters &Hp){
//
//    unsigned int vec;
//    int arg_1, arg_2;
//    int i, ix, iy, iz, ipx, ipy, ipz, nn_i;
//    double E_betanp_temp=0., E_betanm_temp=0.;
//
//    for(iz=0;iz<Lz;iz++){
//        for(iy=0;iy<Ly;iy++){
//            for(ix=0; ix<Lx; ix++){
//                i = ix + Lx * (iy + iz * Ly);
//                for(vec=0; vec<3; vec++) {
//                    ipx=ix;
//                    ipy=iy;
//                    ipz=iz;
//                    if (vec == 0) {
//                        ipx=mod(ix + 1, Lx);
//                    }
//                    if (vec == 1) {
//                        ipy=mod(iy + 1, Ly);
//                    }
//                    if (vec == 2) {
//                        ipz = mod(iz + 1, Lz);
//                    }
//                    nn_i = ipx + Lx * (ipy + ipz * Ly);
//
//                    arg_1 = arg(Site[nn_i].Psi[0] - Site.at(i).Psi[0] - inv_dp*Hp.e * Site[i].A[vec], MaxP);
//                    arg_2 = arg(Site[nn_i].Psi[1] - Site.at(i).Psi[1] - inv_dp*Hp.e * Site[i].A[vec], MaxP);
//                    E_betanm_temp+= vil.potential_bminus.at(OFFSET_POT + arg_1 +MaxP*arg_2);
//                    E_betanp_temp+= vil.potential_bplus.at(OFFSET_POT + arg_1 +MaxP*arg_2);
//                }
//            }
//        }
//    }
//    E_betanp=E_betanp_temp + mis.E_jp + mis.E_gp;
//    E_betanm=E_betanm_temp + mis.E_jp + mis.E_gp;
//}
//
//
//void energy(struct Measures &mis, struct Villain &vil, const std::vector<Node> &Site, struct H_parameters &Hp){
//
//    unsigned int vec;
//    int arg_1, arg_2;
//    int i, ix, iy, iz, nn_i;
//
//    for(iz=0;iz<Lz;iz++){
//        for(iy=0;iy<Ly;iy++){
//            for(ix=0; ix<Lx; ix++){
//                i = ix + Lx * (iy + iz * Ly);
//                for(vec=0; vec<3; vec++) {
//                    if(vec==0){
//                        nn_i= mod(ix+1, Lx) + Lx * (iy + iz * Ly);
//                    }
//                    if(vec==1){
//                        nn_i= ix + Lx * (mod(iy+1,Ly) + iz * Ly);
//                    }
//                    if(vec==2){
//                        nn_i= ix + Lx * (iy + mod(iz+1,Lz) * Ly);
//                    }
//
//                    arg_1 = arg(Site[nn_i].Psi[0] - Site.at(i).Psi[0] - inv_dp*Hp.e * Site[i].A[vec], MaxP);
//                    arg_2 = arg(Site[nn_i].Psi[1] - Site.at(i).Psi[1] - inv_dp*Hp.e * Site[i].A[vec], MaxP);
//                    mis.E+= vil.potential.at(OFFSET_POT + arg_1 +MaxP*arg_2) ;
//                }
//            }
//        }
//    }
//    /****Gauge term***/
//    mis.E+=mis.E_gp;
//    /****Josephson term***/
//    mis.E+=mis.E_jp;
//}

