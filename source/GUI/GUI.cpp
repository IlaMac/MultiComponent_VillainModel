#ifdef GUI_ENABLED

#include "GUI.h"

using namespace std;


GUI::GUI (
  Lattice & lattice
) :
  lattice{lattice}
{
  ////
  //// create window
  ////
  // this->window = new Window{1000, 600, "Villain GUI"};
  this->window = new Window{2*1000, 2*600, "Villain GUI"};


  // ////
  // //// determine length in z-direction and set initial camera position
  // ////
  // if (this->worm.lattice.getNumDimensions() == 1) {

  //   // lattice width
  //   const float latticeWidth = this->worm.lattice.getR(this->worm.lattice.getNumSites() - 1)[0] - this->worm.lattice.getR(0)[0];

  //   // the imaginary time should be two times longer than the width of the lattice
  //   this->z_beg    = 0;
  //   this->z_end    = 2 * latticeWidth;
  //   this->z_length = 2 * latticeWidth;

  //   // determine how far away the camera should be
  //   const float margin = 0.1 * max(this->z_length, latticeWidth);
  //   const float dist_h = 0.5 * (this->z_length + margin)  / (tan( 0.5 * this->window->camera.getHorizontalFOV()));
  //   const float dist_v = 0.5 * (latticeWidth   + margin) / (tan( 0.5 * this->window->camera.getVerticalFOV()));
  //   const auto dist = max(dist_h, dist_v);

  //   // set initial camera position
  //   this->window->camera.setInitialPosition({0.5 * latticeWidth, dist, latticeWidth});
  //   this->window->camera.setInitialFront({0, -1, 0});
  //   this->window->camera.setInitialUp({1, 0, 0});
  //   this->window->camera.reset();

  // } else if (this->worm.lattice.getNumDimensions() == 2) {

  //   // lattice width and height
  //   const float latticeWidth  = this->worm.lattice.getR(this->worm.lattice.getNumSites() - 1)[0] - this->worm.lattice.getR(0)[0];
  //   const float latticeHeight = this->worm.lattice.getR(this->worm.lattice.getNumSites() - 1)[1] - this->worm.lattice.getR(0)[1];

  //   // the imaginary time should be two times longer than max(width, height) of the lattice
  //   this->z_beg    = 0;
  //   this->z_end    = 4 * max(latticeWidth, latticeHeight);
  //   this->z_length = 4 * max(latticeWidth, latticeHeight);

  //   const float phi = glm::radians(45.0f);

  //   const float W = latticeWidth;
  //   const float L = this->z_length;
  //   const float Lc = L * cos(phi);
  //   const float Ls = L * sin(phi);
  //   const float Wc = W * cos(phi);
  //   const float Ws = W * sin(phi);
  //   const float alpha = 0.5 * this->window->camera.getHorizontalFOV();

  //   const float d_h = (Ls + Wc - alpha * (Ws + Lc)) / (2 * alpha);
  //   const float d_v = 0.5 * W * 1.1  / (tan( 0.5 * this->window->camera.getVerticalFOV()));
  //   const float d = max(d_v, d_h);

  //   const float x = (Ls * (d + Ws) - Wc * (d + Lc)) / (2 * d + Lc + Ws);

  //   const glm::vec3 up = {0, 1, 0};
  //   const glm::vec3 front = {sin(phi), 0, cos(phi)};
  //   const glm::vec3 side = glm::cross(front, up);

  //   this->window->camera.setInitialPosition((-d * front + x * side + 0.5f * latticeHeight * up));
  //   this->window->camera.setInitialFront(front);
  //   this->window->camera.setInitialUp({0, 1, 0});
  //   this->window->camera.reset();


  // } else {
  //   cout << "GUI::GUI: cannot handle other than 1D and 2D lattices." << endl;
  // }



  // ////
  // //// set and initialize
  // ////
  // this->worm.displayingWorm = true;
  // this->worm.P_distr.fill(0);




  this->window->camera.setInitialPosition({(Lx - 1)/2.f, (Ly - 1)/2.f, -(Lz + 2)});
  this->window->camera.setInitialFront({0, 0, 1});
  this->window->camera.setInitialUp({0, 1, 0});
  this->window->camera.reset();






  // ////
  // //// listen to key press and key repeat
  // ////
  // window->listenOnKeyPress([] (void* self, const int key, const int scancode, const int mods) {
  //   reinterpret_cast<GUI*>(self)->onKeyPress(key, scancode, mods);
  // }, this);
  // window->listenOnKeyRepeat([] (void* self, const int key, const int scancode, const int mods) {
  //   reinterpret_cast<GUI*>(self)->onKeyRepeat(key, scancode, mods);
  // }, this);


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