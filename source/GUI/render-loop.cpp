#ifdef GUI_ENABLED

#include "GUI.h"

using namespace std;

bool GUI::renderLoop (
  const double & t
) {

  ////
  //// compute velocities from phases
  ////
  this->computeVelocities();


  ////
  //// outline
  ////
  this->window->drawRectangle({0, 0, 0}, {Lx, Ly, Lz}, {0.6, 0.6, 0.6}, true);


  ////
  //// draw sites
  ////
  float r0 = 0.05;
  float r1 = 0.1;
  float l  = 1;


  if (this->actComp < NC) {
    ////
    //// display only a single component
    ////

    for (int iz = 0; iz < Lz; iz++) {
      for (int iy = 0; iy < Ly; iy++) {
        for (int ix = 0; ix < Lx; ix++) {

          // filter correct cross section
          if (this->displayCrossSection) {
            if ( (int) this->crossSectionIndex >= 0       && (int) this->crossSectionIndex < Lx           && ix != (int) this->crossSectionIndex - 0      ) continue;
            if ( (int) this->crossSectionIndex >= Lx      && (int) this->crossSectionIndex < Lx + Ly      && iy != (int) this->crossSectionIndex - Lx     ) continue;
            if ( (int) this->crossSectionIndex >= Lx + Ly && (int) this->crossSectionIndex < Lx + Ly + Lz && iz != (int) this->crossSectionIndex - Lx - Ly) continue;
          }

          auto i = is2i(this->actComp, ix, iy, iz);
          
          const auto & node = this->lattice[is2i(ix, iy, iz)];


          glm::vec3 color = this->colorMapHSV(node.Psi[this->actComp] / (float MaxP));

          const glm::vec3 mid = {ix + 0.5, iy + 0.5, iz + 0.5};

          auto v = this->vs[i] / sqrt(this->maxVs2[this->actComp]);

          const glm::vec3 p0 = mid - 0.4f*v;
          const glm::vec3 p1 = mid + 0.4f*v;

          // this->window->drawSphere(mid, 0.1, color, false);
          this->window->drawArrow(p0, p1, r0, r1, l, color);
        }
      }
    }

  } else {
    ////
    //// display all components
    ////
    float maxV2 = 0;
    for (int a = 0; a < NC; a++) maxV2 = max(this->maxVs2[a], maxV2);

    for (int a = 0; a < NC; a++) {

      glm::vec3 color = this->colorMapHSV(a / (float) NC);

      for (int iz = 0; iz < Lz; iz++) {
        for (int iy = 0; iy < Ly; iy++) {
          for (int ix = 0; ix < Lx; ix++) {

            // filter correct cross section
            if (this->displayCrossSection) {
              if ( (int) this->crossSectionIndex >= 0       && (int) this->crossSectionIndex < Lx           && ix != (int) this->crossSectionIndex - 0      ) continue;
              if ( (int) this->crossSectionIndex >= Lx      && (int) this->crossSectionIndex < Lx + Ly      && iy != (int) this->crossSectionIndex - Lx     ) continue;
              if ( (int) this->crossSectionIndex >= Lx + Ly && (int) this->crossSectionIndex < Lx + Ly + Lz && iz != (int) this->crossSectionIndex - Lx - Ly) continue;
            }

            auto i = is2i(a, ix, iy, iz);

            const glm::vec3 mid = {ix + 0.5, iy + 0.5, iz + 0.5};

            auto v = this->vs[i] / sqrt(maxV2);

            const glm::vec3 p0 = mid - 0.4f*v;
            const glm::vec3 p1 = mid + 0.4f*v;

            // this->window->drawSphere(mid, 0.1, color, false);
            this->window->drawArrow(p0, p1, r0, r1, l, color);
          }
        }
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

  ////
  //// active component
  ////
  {
    x += 200;
    this->window->_drawText("comp: ", x, y, 1, glm::vec3(0.5, 0.5, 0.5));
    const auto str = this->actComp == NC ? "all" : to_string(this->actComp);
    this->window->_drawText(str, x, y, 1, glm::vec3(1, 1, 1));
  }

  ////
  //// beta
  ////
  {
    x += 200;
    std::stringstream stream;
    stream << std::fixed << std::setprecision(5) << this->beta;

    this->window->_drawText("beta: ", x, y, 1, glm::vec3(0.5, 0.5, 0.5));
    this->window->_drawText(stream.str(), x, y, 1, glm::vec3(1, 1, 1));
  }


  return true;
}

#endif