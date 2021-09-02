#include "../Window.h"

using namespace std;


void Window::listenOnKeyRepeat (
  void (* callback) (void *, const int , const int, const int),
  void * self
) {
  this->onKeyRepeatCallbacks.push_back({callback, self});
}