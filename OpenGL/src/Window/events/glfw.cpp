#include "../Window.h"

using namespace std;


////
//// glfs versions
////

void Window::onKey (
  GLFWwindow* window,
  const int key,
  const int scancode,
  const int action,
  const int mods
) {
  auto that = reinterpret_cast<Window*>(glfwGetWindowUserPointer(window));
  if      (action == GLFW_PRESS)   that->onKeyPress(key, scancode, mods);
  else if (action == GLFW_RELEASE) that->onKeyRelease(key, scancode, mods);
  else                             that->onKeyRepeat(key, scancode, mods);
}

void Window::onMouseMove (
  GLFWwindow* window,
  const double xpos,
  const double ypos
) {
  auto that = reinterpret_cast<Window*>(glfwGetWindowUserPointer(window));
  that->onMouseMove(xpos, ypos);
}

void Window::onMouseButton (
  GLFWwindow * window,
  const int button,
  const int action,
  const int mods
) {
  // getting cursor position
  double xpos, ypos;
  glfwGetCursorPos(window, &xpos, &ypos);

  auto that = reinterpret_cast<Window*>(glfwGetWindowUserPointer(window));
  if      (action == GLFW_PRESS)   that->onMousePress(xpos, ypos, button, mods);
  else if (action == GLFW_RELEASE) that->onMouseRelease(xpos, ypos, button, mods);
}

void Window::onScroll (
  GLFWwindow * window,
  const double xoffset,
  const double yoffset
) {
  auto that = reinterpret_cast<Window*>(glfwGetWindowUserPointer(window));
  that->onScroll(xoffset, yoffset);
}