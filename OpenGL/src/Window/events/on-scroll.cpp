#include "../Window.h"

using namespace std;


void Window::onScroll (
  const double & xoffset,
  const double & yoffset
) {
  ////
  //// roll camera
  ////
  double xpos, ypos;
  glfwGetCursorPos(this->glfwWindow, &xpos, &ypos);
  int width, height;
  glfwGetWindowSize(this->glfwWindow, &width, &height);
  this->camera.roll(yoffset, xpos, ypos, width, height);
}