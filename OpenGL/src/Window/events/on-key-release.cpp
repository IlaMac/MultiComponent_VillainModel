#include "../Window.h"

using namespace std;


void Window::onKeyRelease (
  const int key,
  const int scancode,
  const int mods
) {

  ////
  //// camera
  ////
  // deactivate motion
  if      (key == GLFW_KEY_W) this->camera.deactivateForward();
  else if (key == GLFW_KEY_S) this->camera.deactivateBackward();
  else if (key == GLFW_KEY_A) this->camera.deactivateLeftwards();
  else if (key == GLFW_KEY_D) this->camera.deactivateRightwards();
  else if (key == GLFW_KEY_SPACE) this->camera.deactivateUpward();
  else if (key == GLFW_KEY_LEFT_CONTROL) this->camera.deactivateDownward();
  else if (key == GLFW_KEY_LEFT_SHIFT) this->camera.deactivateBoost();

}