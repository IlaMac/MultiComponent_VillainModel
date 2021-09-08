#ifdef GUI_ENABLED

#include "GUI.h"
#include <glm/gtx/norm.hpp>

using namespace std;


static inline
float
grad (
  int theta1,
  int theta0
) {
  // average difference of two neighboring bonds
  return (theta1 - theta0) * M_PI / MaxP;
}


void GUI::computeVelocities () {

  for (int a = 0; a < NC; a++) {
    this->maxVs2[a] = 0;

    for (int iz = 0; iz < Lz; iz++) {
      for (int iy = 0; iy < Ly; iy++) {
        for (int ix = 0; ix < Lx; ix++) {
          const auto i = is2i(a, ix, iy, iz);

          const auto ix0 = is2i(ix == 0      ? Lx - 1 : ix - 1, iy, iz);
          const auto ix1 = is2i(ix == Lx - 1 ? 0      : ix + 1, iy, iz);
          this->vs[i].x = grad(this->lattice[ix1].Psi[a], this->lattice[ix0].Psi[a]);

          const auto iy0 = is2i(ix, iy == 0      ? Ly - 1 : iy - 1, iz);
          const auto iy1 = is2i(ix, iy == Ly - 1 ? 0      : iy + 1, iz);
          this->vs[i].y = grad(this->lattice[iy1].Psi[a], this->lattice[iy0].Psi[a]);

          const auto iz0 = is2i(ix, iy, iz == 0      ? Lz - 1 : iz - 1);
          const auto iz1 = is2i(ix, iy, iz == Lz - 1 ? 0      : iz + 1);
          this->vs[i].z = grad(this->lattice[iz1].Psi[a], this->lattice[iz0].Psi[a]);

          this->maxVs2[a] = max(glm::length2(this->vs[i]), this->maxVs2[a]);
        }
      }
    }
  }
}


#endif