//
// Created by ilaria on 2020-12-15.
//

#include "pt.h"
#include "mpi_tools.h"
#include <mpi.h>

void initialize_PTarrays(struct PT_parameters &PTp, struct PTroot_parameters &PTroot, struct H_parameters &Hp) {
    int p;
    double beta_low, beta_high, delta_beta;

    if (Hp.b_high > Hp.b_low) { //Paranoic check
        beta_low = Hp.b_low;
        beta_high = Hp.b_high;
    } else {
        beta_low = Hp.b_high;
        beta_high = Hp.b_low;
    }
    PTroot.beta.resize(PTp.np, 0.0);
    PTroot.beta_p.resize(PTp.np, 0.0);
    PTroot.beta_m.resize(PTp.np, 0.0);
    //PTroot.All_Energies.resize(PTp.np, 0.0);
    PTroot.E_rank_beta.resize(PTp.np, 0.0);
    PTroot.E_rank_betap.resize(PTp.np, 0.0);
    PTroot.E_rank_betam.resize(PTp.np, 0.0);

    PTroot.Villain_beta.resize(PTp.np);

    PTroot.ind_to_rank.resize(PTp.np, 0);
    PTroot.rank_to_ind.resize(PTp.np, 0);
    delta_beta = (beta_high - beta_low) / (PTp.np - 1);
    for (p = 0; p < PTp.np; p++) {
        PTroot.rank_to_ind[p] = p;
        PTroot.ind_to_rank[p] = p;
        PTroot.beta[p] = beta_low + p * delta_beta;
    }
    for (p = 0; p < PTp.np; p++) {
        PTroot.beta_p[p] = PTroot.beta[ (PTp.np+p+1)%PTp.np ];
        PTroot.beta_m[p] = PTroot.beta[ (PTp.np+p-1)%PTp.np ];
    }
}
void parallel_temp(double &my_E , double &E_betanp, double &E_betanm, double &beta_p, double &beta_m,  double &my_beta, int &my_ind, struct Villain &vil, struct PT_parameters &PTp, struct PTroot_parameters &PTroot){

    double coin;
    double n_rand;
    //, delta_E, delta_beta;
    double Delta;
    double E_rank_betann, E_ranknn_beta;
    double oldbeta_i, oldbeta_nn;

    int i=0,p, nn=0, ind_nn;
    int oldrank_i, oldrank_nn;
    int newrank_i, newrank_nn;

    struct Villain vill_oldbeta;
    struct Villain vill_oldbeta_nn;

    std::vector<double> temp_potential(vil.potential.size()*PTp.np);
    //Create a datatype for the nth worker_results[n] struct
    MPI_Datatype vill_struct;
    MPI_Type_contiguous(9*MaxP*MaxP,MPI_DOUBLE,&vill_struct);
    MPI_Type_commit(&vill_struct);

    MPI_Gather(&my_E, 1, MPI_DOUBLE, PTroot.E_rank_beta.data(), 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Gather(&E_betanp, 1, MPI_DOUBLE, PTroot.E_rank_betap.data(), 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Gather(&E_betanm, 1, MPI_DOUBLE, PTroot.E_rank_betam.data(), 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Gather(&beta_m, 1, MPI_DOUBLE, PTroot.beta_m.data(), 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Gather(&beta_p, 1, MPI_DOUBLE, PTroot.beta_p.data(), 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
//    MPI_Gather(&vil, 1 , vill_struct, PTroot.Villain_beta.data(), 1, vill_struct, PTp.root, MPI_COMM_WORLD);

    auto tmp_potential        = mpi::gather(vil.potential         ,PTp.rank, PTp.root,PTp.np, MPI_DOUBLE);
    auto tmp_potential_bplus  = mpi::gather(vil.potential_bplus   ,PTp.rank, PTp.root,PTp.np, MPI_DOUBLE);
    auto tmp_potential_bminus = mpi::gather(vil.potential_bminus  ,PTp.rank, PTp.root,PTp.np, MPI_DOUBLE);
    auto tmp_upotential       = mpi::gather(vil.upotential        ,PTp.rank, PTp.root,PTp.np, MPI_DOUBLE);
    auto tmp_d1_potential     = mpi::gather(vil.d1_potential      ,PTp.rank, PTp.root,PTp.np, MPI_DOUBLE);
    auto tmp_d2_potential     = mpi::gather(vil.d2_potential      ,PTp.rank, PTp.root,PTp.np, MPI_DOUBLE);
    auto tmp_d11_potential    = mpi::gather(vil.d11_potential     ,PTp.rank, PTp.root,PTp.np, MPI_DOUBLE);
    auto tmp_d22_potential    = mpi::gather(vil.d22_potential     ,PTp.rank, PTp.root,PTp.np, MPI_DOUBLE);
    auto tmp_d12_potential    = mpi::gather(vil.d12_potential     ,PTp.rank, PTp.root,PTp.np, MPI_DOUBLE);
    if(PTp.rank== PTp.root){
        for(size_t b = 0; b < PTroot.Villain_beta.size(); b++){
            PTroot.Villain_beta[b].potential = tmp_potential[b];
            PTroot.Villain_beta[b].potential_bplus = tmp_potential_bplus[b];
            PTroot.Villain_beta[b].potential_bminus = tmp_potential_bminus[b];
            PTroot.Villain_beta[b].upotential = tmp_upotential[b];
            PTroot.Villain_beta[b].d1_potential = tmp_d1_potential[b];
            PTroot.Villain_beta[b].d2_potential = tmp_d2_potential[b];
            PTroot.Villain_beta[b].d11_potential = tmp_d11_potential[b];
            PTroot.Villain_beta[b].d22_potential = tmp_d22_potential[b];
            PTroot.Villain_beta[b].d12_potential = tmp_d12_potential[b];
        }
    }


    if (PTp.rank == PTp.root) { //Root forms the pairs and decides (given the energies and the betas) which pairs will swap
        //Pair Formation
        coin = rn::uniform_real_box(0, 1);
        if (coin < 0.5) { //each even rank wil be paired with its right neighbour
            nn = +1;
        } else if (coin >= 0.5) { //each even rank wil be paired with its left neighbour
            nn = -1;
        }
        while (i < PTp.np) {
            ind_nn = (PTp.np + i + nn) % PTp.np;
            oldrank_i = PTroot.ind_to_rank[i];
            oldrank_nn = PTroot.ind_to_rank[ind_nn];
            if (nn == +1) {
                E_rank_betann = PTroot.E_rank_betap[oldrank_i]; /*H[x_i, \beta_nn]*/
                E_ranknn_beta = PTroot.E_rank_betam[oldrank_nn]; /*H[x_nn, \beta_i]*/
            } else if (nn == -1) {
                E_rank_betann = PTroot.E_rank_betam[oldrank_i]; /*H[x_i, \beta_nn]*/
                E_ranknn_beta = PTroot.E_rank_betap[oldrank_nn]; /*H[x_nn, \beta_i]*/
            }

            /*Delta= beta_nn*[H(x_nn, beta_nn) - H(x_i, beta_nn)] - beta_i*[H(x_nn,beta_i) - H(x_i, beta_i)] */
            Delta = PTroot.beta[oldrank_nn] * (PTroot.E_rank_beta[oldrank_nn]- E_rank_betann) - PTroot.beta[oldrank_i] * (E_ranknn_beta  - PTroot.E_rank_beta[oldrank_i]);
            //swapping condition
            n_rand = rn::uniform_real_box(0, 1);
            if ((n_rand < exp(Delta)) ) {

                //swap indices in the rank_to_ind array
                PTroot.rank_to_ind[oldrank_i] = ind_nn;
                PTroot.rank_to_ind[oldrank_nn] = i;

                //swap indices in the ind_to_rank array
                newrank_i = oldrank_nn;
                PTroot.ind_to_rank[i] = newrank_i;
                newrank_nn = oldrank_i;
                PTroot.ind_to_rank[ind_nn] = newrank_nn;
                //swap beta
                oldbeta_i = PTroot.beta[oldrank_i];
                oldbeta_nn = PTroot.beta[oldrank_nn];
                PTroot.beta[oldrank_i] = oldbeta_nn;
                PTroot.beta[oldrank_nn] = oldbeta_i;
                //swap vill struct
                vill_oldbeta=PTroot.Villain_beta[oldrank_i];
                vill_oldbeta_nn=PTroot.Villain_beta[oldrank_nn];
                PTroot.Villain_beta[oldrank_i]=vill_oldbeta_nn;
                PTroot.Villain_beta[oldrank_nn]=vill_oldbeta;
            }
            i += 2;
        }

        for (p = 0; p < PTp.np; p++) {
            PTroot.beta_p[PTroot.ind_to_rank[p]] = PTroot.beta[ PTroot.ind_to_rank[(PTp.np+p+1)%PTp.np] ];
            PTroot.beta_m[PTroot.ind_to_rank[p]] = PTroot.beta[ PTroot.ind_to_rank[(PTp.np+p-1)%PTp.np ]];
        }

    }
    MPI_Scatter(PTroot.beta.data(), 1, MPI_DOUBLE, &my_beta, 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Scatter(PTroot.rank_to_ind.data(), 1, MPI_INT, &my_ind, 1, MPI_INT, PTp.root, MPI_COMM_WORLD);
    MPI_Scatter(PTroot.beta_m.data(), 1, MPI_DOUBLE, &beta_m, 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Scatter(PTroot.beta_p.data(), 1, MPI_DOUBLE, &beta_p, 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
   // MPI_Scatter(PTroot.Villain_beta.data(), 1, vill_struct, &vil, 1, vill_struct, PTp.root, MPI_COMM_WORLD);

    if(PTp.rank== PTp.root){
        for(size_t i = 0; i < PTroot.Villain_beta.size(); i++){
            tmp_potential[i] = PTroot.Villain_beta[i].potential;
            tmp_potential_bplus[i] = PTroot.Villain_beta[i].potential_bplus;
            tmp_potential_bminus[i] = PTroot.Villain_beta[i].potential_bminus;
            tmp_upotential[i] = PTroot.Villain_beta[i].upotential;
            tmp_d1_potential[i] = PTroot.Villain_beta[i].d1_potential;
            tmp_d2_potential[i] = PTroot.Villain_beta[i].d2_potential;
            tmp_d11_potential[i] = PTroot.Villain_beta[i].d11_potential;
            tmp_d22_potential[i] = PTroot.Villain_beta[i].d22_potential;
            tmp_d12_potential[i] = PTroot.Villain_beta[i].d12_potential;
        }
    }
//    vil.potential = mpi::scatter(merged_potential, PTp.root, PTp.np, MPI_DOUBLE );
    mpi::scatter(tmp_potential,vil.potential, PTp.root, PTp.rank, PTp.np, MPI_DOUBLE );
    mpi::scatter(tmp_potential_bplus,vil.potential_bplus, PTp.root, PTp.rank, PTp.np, MPI_DOUBLE );
    mpi::scatter(tmp_potential_bminus,vil.potential_bminus, PTp.root, PTp.rank, PTp.np, MPI_DOUBLE );
    mpi::scatter(tmp_upotential,vil.upotential, PTp.root, PTp.rank, PTp.np, MPI_DOUBLE );
    mpi::scatter(tmp_d1_potential,vil.d1_potential, PTp.root, PTp.rank, PTp.np, MPI_DOUBLE );
    mpi::scatter(tmp_d2_potential,vil.d2_potential, PTp.root, PTp.rank, PTp.np, MPI_DOUBLE );
    mpi::scatter(tmp_d11_potential,vil.d11_potential, PTp.root, PTp.rank, PTp.np, MPI_DOUBLE );
    mpi::scatter(tmp_d22_potential,vil.d22_potential, PTp.root, PTp.rank, PTp.np, MPI_DOUBLE );
    mpi::scatter(tmp_d12_potential,vil.d12_potential, PTp.root, PTp.rank, PTp.np, MPI_DOUBLE );

}
