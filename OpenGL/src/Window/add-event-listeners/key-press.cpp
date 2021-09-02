#include "../Window.h"

using namespace std;


void Window::listenOnKeyPress (
  void (* callback) (void *, const int , const int, const int),
  void * self
) {
  this->onKeyPressCallbacks.push_back({callback, self});
}