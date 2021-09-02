#include "../Window.h"

using namespace std;


void Window::onKeyPress (
  const int key,
  const int scancode,
  const int mods
) {
  ////
  //// close window on escape
  ////
  if (key == GLFW_KEY_ESCAPE) glfwSetWindowShouldClose(this->glfwWindow, GLFW_TRUE);


  ////
  //// camera
  ////

  // activate motion
  if      (key == GLFW_KEY_W) this->camera.activateForward();
  else if (key == GLFW_KEY_S) this->camera.activateBackward();
  else if (key == GLFW_KEY_A) this->camera.activateLeftwards();
  else if (key == GLFW_KEY_D) this->camera.activateRightwards();
  else if (key == GLFW_KEY_SPACE) this->camera.activateUpward();
  else if (key == GLFW_KEY_LEFT_CONTROL) this->camera.activateDownward();
  else if (key == GLFW_KEY_LEFT_SHIFT) this->camera.activateBoost();

  // reset camera
  if (key == GLFW_KEY_BACKSPACE) this->camera.reset();


  ////
  //// execute listener callbacks
  ////
  for (auto [callback, self] : this->onKeyPressCallbacks) {
    callback(self, key, scancode, mods);
  }
}