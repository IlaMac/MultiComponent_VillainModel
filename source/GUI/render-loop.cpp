#ifdef GUI_ENABLED

#include "GUI.h"

#include <glm/gtx/norm.hpp>

using namespace std;
// using std::chrono::system_clock;

bool GUI::renderLoop (
  const double & t
) {



  for (int iz = 0; iz < Lz; iz++) {
    for (int iy = 0; iy < Ly; iy++) {
      for (int ix = 0; ix < Lx; ix++) {
        unsigned i = ix + Lx * (iy + iz * Ly);

        const auto & node = this->lattice[i];
        // cout << node.Psi[0] << endl; // (float MaxP) << endl;
        this->window->drawSphere({ix, iy, iz}, 0.5, {1, 0 ,0}, false);
      }
    }
  }






  // this->window->drawSphere({0, 0, 8}, 1, {0, 1 ,0}, false);




  return true;
}

#endif