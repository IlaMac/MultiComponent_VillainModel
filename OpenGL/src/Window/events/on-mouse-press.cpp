#include "../Window.h"

using namespace std;


void Window::onMousePress (
  const double & xpos,
  const double & ypos,
  const int button,
  const int mods
) {
  ////
  //// orient camera
  ////
  if (button == GLFW_MOUSE_BUTTON_1) {
    // // get window size
    // int width, height;
    // glfwGetWindowSize(this->glfwWindow, &width, &height);

    this->camera.activateOrientation(xpos, ypos);//, width, height);
  }
}