#include "../Window.h"

using namespace std;


void Window::onMouseMove (
  const double & xpos,
  const double & ypos
) {
  ////
  //// orient camera
  ////
  this->camera.orient(xpos, ypos);
}