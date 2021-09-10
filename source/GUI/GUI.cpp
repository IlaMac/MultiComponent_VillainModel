#ifdef GUI_ENABLED

#include "GUI.h"

using namespace std;


GUI::GUI (
  const Lattice        & lattice,
  struct MC_parameters & MCp,
  struct H_parameters  & Hp,
  const double           beta,
  struct Villain       & villain
) :
  lattice{lattice},
  MCp{MCp},
  Hp{Hp},
  beta{beta},
  villain{villain}
{
  ////
  //// create window
  ////
  // this->window = new Window{1000, 600, "Villain GUI"};
  this->window = new Window{2*1000, 2*600, "Villain GUI"};


  ////
  //// listen to key press and key repeat
  ////
  window->listenOnKeyPress([] (void* self, const int key, const int scancode, const int mods) {
    reinterpret_cast<GUI*>(self)->onKeyPress(key, scancode, mods);
  }, this);
  window->listenOnKeyRepeat([] (void* self, const int key, const int scancode, const int mods) {
    reinterpret_cast<GUI*>(self)->onKeyRepeat(key, scancode, mods);
  }, this);


  ////
  //// initiate camera
  ////
  const float latticeWidth  = Lx;
  const float latticeHeight = Ly;
  const float latticeDepth  = Lz;

  const float phi = glm::radians(35.0f);

  const float W = latticeWidth;
  const float L = latticeDepth;
  const float Lc = L * cos(phi);
  const float Ls = L * sin(phi);
  const float Wc = W * cos(phi);
  const float Ws = W * sin(phi);
  const float alpha = 0.5 * this->window->camera.getHorizontalFOV();

  const float d_h = (Ls + Wc - alpha * (Ws + Lc)) / (2 * alpha);
  const float d_v = 0.5 * W * 1.1  / (tan( 0.5 * this->window->camera.getVerticalFOV()));
  const float d = max(d_v, d_h);

  const float x = (Ls * (d + Ws) - Wc * (d + Lc)) / (2 * d + Lc + Ws);

  const glm::vec3 up = {0, 1, 0};
  const glm::vec3 front = {sin(phi), 0, cos(phi)};
  const glm::vec3 side = glm::cross(front, up);

  this->window->camera.setInitialPosition((-d * front + x * side + 0.5f * latticeHeight * up));
  this->window->camera.setInitialFront(front);
  this->window->camera.setInitialUp({0, 1, 0});
  this->window->camera.reset();


  ////
  //// setup variables
  ////
  this->vs.resize(Lx * Ly * Lz * NC);

  // individual (3) + differences ( N*[N-1]/2 ) + sum
  this->avgVsLength.resize(NC + NC*(NC-1)/2 + 1);


  ////
  //// start main render loop
  ////
  window->show([] (void* self, const double & t) {
    return reinterpret_cast<GUI*>(self)->renderLoop(t);
  }, this);
}



GUI::~GUI () {
  // remove the dynamically allocated Window object
  delete this->window;
}

#endif