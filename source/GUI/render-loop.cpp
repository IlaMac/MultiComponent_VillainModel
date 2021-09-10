#ifdef GUI_ENABLED

#include "GUI.h"

using namespace std;

#define DRAW_ARROW(mid, v, r0, r1, l, color)                             \
{                                                                        \
  this->window->drawArrow(mid - 0.4f*v, mid + 0.4f*v, r0, r1, l, color); \
}
// this->window->drawSphere(mid, 0.1, color, false);

bool GUI::renderLoop (
  const double & t
) {
  // cout << "sum + difference + print avg velocity" << endl;

  ////
  //// compute velocities from phases
  ////
  this->computeVelocities();


  ////
  //// outline
  ////
  this->window->drawRectangle({0, 0, 0}, {Lx, Ly, Lz}, {0.6, 0.6, 0.6}, true);


  ////
  //// draw velocities
  ////
  //// TODO: maybe it is better to compute all sum and differences of velocities in computeVelocities()? rather than here
  ////
  float r0 = 0.05;
  float r1 = 0.1;
  float l  = 1;



  for (int iz = 0; iz < Lz; iz++) {
    for (int iy = 0; iy < Ly; iy++) {
      for (int ix = 0; ix < Lx; ix++) {

        // filter correct cross section
        if (this->displayCrossSection) {
          if ( (int) this->crossSectionIndex >= 0       && (int) this->crossSectionIndex < Lx           && ix != (int) this->crossSectionIndex - 0      ) continue;
          if ( (int) this->crossSectionIndex >= Lx      && (int) this->crossSectionIndex < Lx + Ly      && iy != (int) this->crossSectionIndex - Lx     ) continue;
          if ( (int) this->crossSectionIndex >= Lx + Ly && (int) this->crossSectionIndex < Lx + Ly + Lz && iz != (int) this->crossSectionIndex - Lx - Ly) continue;
        }

        // reference to node
        const auto & node = this->lattice[is2i(ix, iy, iz)];

        // middle point of lattice
        const glm::vec3 mid = {ix + 0.5, iy + 0.5, iz + 0.5};

        glm::vec3 v;
        glm::vec3 color;

        if (this->actComp < NC) {
          // individual components
          auto i = is2i(this->actComp, ix, iy, iz);

          v     = this->vs[i] / this->avgVsLength[this->actComp];
          color = this->colorMapHSV(node.Psi[this->actComp] / (float MaxP));

          DRAW_ARROW(mid, v, r0, r1, l, color);
        }
        else if (this->actComp == NC) {
          // superimposed components
          for (unsigned a = 0; a < NC; a++) {
            auto i = is2i(a, ix, iy, iz);

            v     = this->vs[i] / this->avgVsLength[a];
            color = this->colorMapHSV(a / (float) NC);

            DRAW_ARROW(mid, v, r0, r1, l, color);
          }
        }
        else if (this->actComp < NC + 1 + NC*(NC-1)/2) {
          // component differences
          const auto [a, b] = i2iaib(this->actComp - NC - 1);

          auto i = is2i(a, ix, iy, iz);
          auto j = is2i(b, ix, iy, iz);

          v     = (this->vs[i] - this->vs[j]) / this->avgVsLength[this->actComp - NC + 1];
          color = this->colorMapHSV(arg(node.Psi[a] - node.Psi[b], MaxP) / (float MaxP));

          DRAW_ARROW(mid, v, r0, r1, l, color);
        }
        else {
          // component sum
          glm::vec3 sum_v = {0, 0, 0};
          int       sum_p = 0;

          // superimposed components
          for (unsigned a = 0; a < NC; a++) {
            auto i = is2i(a, ix, iy, iz);

            sum_v += this->vs[i];
            sum_p += node.Psi[a];
          }

          sum_v /= this->avgVsLength.back();
          color = this->colorMapHSV(arg(sum_p, MaxP) / (float MaxP));

          DRAW_ARROW(mid, sum_v, r0, r1, l, color);
        }
      }
    }
  }




  ////
  //// display average velocity length
  ////
  {
    float x = 10;
    float y = this->window->getHeight() - 40;


    this->window->_drawText("<|v|>: ", x, y, 1, glm::vec3(0.5, 0.5, 0.5));


    std::stringstream stream;
    stream << std::fixed << std::setprecision(3);
    if (this->actComp < NC) {
      // individual components
      stream << this->avgVsLength[this->actComp];
      this->window->_drawText(stream.str(), x, y, 1, glm::vec3(1, 1, 1));
    }
    else if (this->actComp == NC) {
      // superimposed components
      for (unsigned a = 0; a < NC; a++) {
        stream.str("");
        stream << this->avgVsLength[a];
        glm::vec3 color = this->colorMapHSV(a / (float) NC);
        this->window->_drawText(stream.str(), x, y, 1, color);
        if (a < NC - 1)
          this->window->_drawText(", ", x, y, 1, glm::vec3(0.5, 0.5, 0.5));
      }
    }
    else if (this->actComp < NC + 1 + NC*(NC-1)/2) {
      // component differences
      stream << this->avgVsLength[this->actComp - NC + 1];
      this->window->_drawText(stream.str(), x, y, 1, glm::vec3(1, 1, 1));
    }
    else {
      // component sum
      stream << this->avgVsLength.back();
      this->window->_drawText(stream.str(), x, y, 1, glm::vec3(1, 1, 1));
    }

  }


  // cout << this->avgVsLength << endl;



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
  //// beta
  ////
  {
    x += 200;
    std::stringstream stream;
    stream << std::fixed << std::setprecision(5) << this->beta;

    this->window->_drawText("beta: ", x, y, 1, glm::vec3(0.5, 0.5, 0.5));
    this->window->_drawText(stream.str(), x, y, 1, glm::vec3(1, 1, 1));
  }

  ////
  //// active component
  ////
  {
    x += 200;
    this->window->_drawText("comp: ", x, y, 1, glm::vec3(0.5, 0.5, 0.5));
    float _x = x, _y = y;


    if (this->actComp < NC) {
      // individual components
      const auto str = to_string(this->actComp);
      this->window->_drawText(str, _x, _y, 1, glm::vec3(1, 1, 1));
    }
    else if (this->actComp == NC) {
      // superimposed components
      for (unsigned a = 0; a < NC; a++) {
        glm::vec3 color = this->colorMapHSV(a / (float) NC);
        this->window->_drawText(to_string(a), _x, _y, 1, color);
        if (a < NC - 1)
          this->window->_drawText(", ", _x, _y, 1, glm::vec3(0.5, 0.5, 0.5));
      }
    }
    else if (this->actComp < NC + 1 + NC*(NC-1)/2) {
      // component differences
      const auto [a, b] = i2iaib(this->actComp - NC - 1);
      const auto str = to_string(a) + " - " + to_string(b);
      this->window->_drawText(str, _x, _y, 1, glm::vec3(1, 1, 1));
    }
    else {
      // component sum
      string str;
      for (unsigned a = 0; a < NC; a++) {
        str += to_string(a);
        if (a < NC - 1) str += " + ";
      }
      this->window->_drawText(str, _x, _y, 1, glm::vec3(1, 1, 1));
    }

  }


  return true;
}

#endif