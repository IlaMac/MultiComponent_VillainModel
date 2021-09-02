#include "../Window.h"

using namespace std;


void Window::onMouseRelease (
  const double & xpos,
  const double & ypos,
  const int button,
  const int mods
) {
  ////
  //// orient camera
  ////
  if (button == GLFW_MOUSE_BUTTON_1) this->camera.deactivateOrientation();
}