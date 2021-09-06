#ifdef GUI_ENABLED

#include "GUI.h"

#include <glm/gtx/norm.hpp>

using namespace std;
// using std::chrono::system_clock;

bool GUI::renderLoop (
  const double & t
) {


  ////
  //// outline
  ////
  this->window->drawRectangle({-0.5, -0.5, -0.5}, {Lx, Ly, Lz}, {0.6, 0.6, 0.6}, true);



  ////
  //// draw sites
  ////
  float r0 = 0.05;
  float r1 = 0.1;
  float l  = 0.15;

  for (int iz = 0; iz < Lz; iz++) {
    for (int iy = 0; iy < Ly; iy++) {
      for (int ix = 0; ix < Lx; ix++) {

        // filter correct cross section
        if (this->displayCrossSection) {
          if ( (int) this->crossSectionIndex >= 0       && (int) this->crossSectionIndex < Lx           && ix != (int) this->crossSectionIndex - 0      ) continue;
          if ( (int) this->crossSectionIndex >= Lx      && (int) this->crossSectionIndex < Lx + Ly      && iy != (int) this->crossSectionIndex - Lx     ) continue;
          if ( (int) this->crossSectionIndex >= Lx + Ly && (int) this->crossSectionIndex < Lx + Ly + Lz && iz != (int) this->crossSectionIndex - Lx - Ly) continue;
        }

        // unsigned i = ix + Lx * (iy + iz * Ly);
        // const auto & node = this->lattice[i];
        // cout << node.Psi[0] << endl; // (float MaxP) << endl;


        glm::vec3 color = {0.5 + 0.5*ix/(Lx-1.), 0.5 + 0.5*iy/(Ly-1.), 0.5 + 0.5*iz/(Lz-1.)};


        const glm::vec3 mid = {ix, iy, iz};

        const glm::vec3 v = {0, cos(ix/(Lx-1.) * 2*M_PI), sin(ix/(Lx-1.) * 2*M_PI)};

        const glm::vec3 p0 = mid - 0.4f*v;
        const glm::vec3 p1 = mid + 0.4f*v;

        // this->window->drawSphere(mid, 0.1, color, false);
        this->window->drawArrow(p0, p1, r0, r1, l, color);
      }
    }
  }









  ////
  //// displayframerate
  ////
  float x = 10;
  float y = 10;
  {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(0) << this->window->getFramerate();

    this->window->_drawText("fps: ", x, y, 1, glm::vec3(0.5, 0.5, 0.5));
    float _x = x, _y = y;
    this->window->_drawText(stream.str(), _x, _y, 1, glm::vec3(1, 1, 1));
  }


  return true;
}

#endif