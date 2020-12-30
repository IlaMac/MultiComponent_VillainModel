//
// Created by ilaria on 2020-12-15.
//

#include "pt.h"

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
void parallel_temp(double &my_E , double &E_betanp, double &E_betanm,  double &my_beta, int &my_ind, struct PT_parameters &PTp, struct PTroot_parameters &PTroot){

    double coin;
    double n_rand;
    //, delta_E, delta_beta;
    double Delta;
    double E_rank_betann, E_ranknn_beta;
    double oldbeta_i, oldbeta_nn;
    int i=0, nn=0, ind_nn;
    int oldrank_i, oldrank_nn;
    int newrank_i, newrank_nn;


    MPI_Gather(&my_E, 1, MPI_DOUBLE, PTroot.E_rank_beta.data(), 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Gather(&E_betanp, 1, MPI_DOUBLE, PTroot.E_rank_betap.data(), 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Gather(&E_betanm, 1, MPI_DOUBLE, PTroot.E_rank_betam.data(), 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);

    if (PTp.rank == PTp.root) { //Root forms the pairs and decides (given the energies and the betas) which pairs will swap
        //Pair Formation
        coin = rn::uniform_real_box(0, 1);
        if (coin < 0.5) { //each even rank wil be paired with its right neighbour
            nn = +1;
        } else if (coin >= 0.5) { //each even rank wil be paired with its left neighbour
            nn = -1;
        }
        while (i < PTp.np) {
            n_rand = rn::uniform_real_box(0, 1);
            ind_nn = (PTp.np + i + nn) % PTp.np;
            oldrank_i = PTroot.ind_to_rank[i];
            oldrank_nn = PTroot.ind_to_rank[ind_nn];
            if (nn == 1) {
                E_rank_betann = PTroot.E_rank_betap[oldrank_i];
                E_ranknn_beta = PTroot.E_rank_betam[oldrank_nn];
            } else if (nn == -1) {
                E_rank_betann = PTroot.E_rank_betam[oldrank_i];
                E_ranknn_beta = PTroot.E_rank_betap[oldrank_nn];
            }

            /*Delta= beta_new*[H(x, beta_new) - H(x_new, beta_new)] - beta*[H(x,beta) - H(x_new, beta)] */
            Delta = PTroot.beta[oldrank_nn] * (E_rank_betann - PTroot.E_rank_beta[oldrank_nn]) -
                    PTroot.beta[oldrank_i] * (PTroot.E_rank_beta[oldrank_i] - E_ranknn_beta);

            //swapping condition
            if (n_rand < exp(-Delta)) {

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
            }
            i += 2;
        }
        for (i = 0; i < PTp.np; i++) {
            PTroot.beta_p[i] = PTroot.beta[ (PTp.np+i+1)%PTp.np ];
            PTroot.beta_m[i] = PTroot.beta[ (PTp.np+i-1)%PTp.np ];
        }
    }
    MPI_Scatter(PTroot.beta.data(), 1, MPI_DOUBLE, &my_beta, 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Scatter(PTroot.rank_to_ind.data(), 1, MPI_INT, &my_ind, 1, MPI_INT, PTp.root, MPI_COMM_WORLD);

}