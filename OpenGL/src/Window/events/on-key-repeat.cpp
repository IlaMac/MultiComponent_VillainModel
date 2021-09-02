#include "../Window.h"

using namespace std;


void Window::onKeyRepeat(
  const int key,
  const int scancode,
  const int mods
) {
  ////
  //// execute listener callbacks
  ////
  for (auto [callback, self] : this->onKeyRepeatCallbacks) {
    callback(self, key, scancode, mods);
  }
}