#ifdef GUI_ENABLED

#include "GUI.h"
#include <glm/gtx/norm.hpp>

using namespace std;


static inline
float
grad (
  int theta1,
  int theta,
  int theta0
) {
  // at most the phases point in opposite directions, i.e., a maximum difference of ± 180 rather than ± 360
  int diff1 = arg(theta1 - theta,  MaxP);
  int diff0 = arg(theta  - theta0, MaxP);

  // average difference of two neighboring bonds
  return 0.5*(diff1 + diff0) * 2*M_PI/MaxP;
}


void GUI::computeVelocities () {
  ////
  //// TODO: compute only if not already done! store some sort of update number count which is incremented when when an update is carried out
  ////


  ////
  //// compute velocities
  ////
  for (int a = 0; a < NC; a++) {
    for (int iz = 0; iz < Lz; iz++) {
      for (int iy = 0; iy < Ly; iy++) {
        for (int ix = 0; ix < Lx; ix++) {
          const auto i = is2i(a, ix, iy, iz);

          const auto ix0 = is2i(ix == 0      ? Lx - 1 : ix - 1, iy, iz);
          const auto ix1 = is2i(ix == Lx - 1 ? 0      : ix + 1, iy, iz);
          this->vs[i].x = grad(this->lattice[ix1].Psi[a], this->lattice[ix].Psi[a], this->lattice[ix0].Psi[a]);

          const auto iy0 = is2i(ix, iy == 0      ? Ly - 1 : iy - 1, iz);
          const auto iy1 = is2i(ix, iy == Ly - 1 ? 0      : iy + 1, iz);
          this->vs[i].y = grad(this->lattice[iy1].Psi[a], this->lattice[iy].Psi[a], this->lattice[iy0].Psi[a]);

          const auto iz0 = is2i(ix, iy, iz == 0      ? Lz - 1 : iz - 1);
          const auto iz1 = is2i(ix, iy, iz == Lz - 1 ? 0      : iz + 1);
          this->vs[i].z = grad(this->lattice[iz1].Psi[a], this->lattice[iz].Psi[a], this->lattice[iz0].Psi[a]);
        }
      }
    }
  }


  ////
  //// compute averages
  ////

  // reset
  fill(this->avgVsLength.begin(), this->avgVsLength.end(), 0);

  for (int iz = 0; iz < Lz; iz++) {
    for (int iy = 0; iy < Ly; iy++) {
      for (int ix = 0; ix < Lx; ix++) {

        glm::vec3 sum = {0, 0, 0};

        for (int a = 0; a < NC; a++) {
          const auto i = is2i(a, ix, iy, iz);

          // add
          sum += this->vs[i];

          // individual
          this->avgVsLength[a] += glm::length(this->vs[i]);

          // difference
          for (int b = a + 1; b < NC; b++) {
            const auto j = is2i(b, ix, iy, iz);

            const auto I = iaib2i(a, b);
            this->avgVsLength[NC + I] += glm::length(this->vs[i] - this->vs[j]);
          }
        }

        // sum
        this->avgVsLength[NC + NC*(NC-1)/2] += glm::length(sum);
      }
    }
  }

  // normalize
  for (auto & x : this->avgVsLength) x /= Lz * Ly * Lx;
}


#endif