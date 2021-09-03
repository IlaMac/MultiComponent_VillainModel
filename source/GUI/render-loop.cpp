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
        // this->window->drawSphere({ix, iy, iz}, 0.5, {1, 0 ,0}, false);
        this->window->drawRectangle({ix-0.5, iy-0.5, iz-0.5}, {1, 1, 1}, {ix/(Lx-1.), iy/(Ly-1.), iz/(Lz-1.)}, false);
      }
    }
  }





  return true;
}

#endif