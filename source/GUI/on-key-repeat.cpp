#ifdef GUI_ENABLED

#include "GUI.h"

using namespace std;


void GUI::onKeyRepeat (
  const int key,
  const int scancode,
  const int mods
) {
  ////
  //// update configuration
  ////
  if (key == GLFW_KEY_ENTER) {
    // press enter -> single updates
    this->singleUpdate();

    // inactivate multiple updates process
    // (the thread listens to "this->multipleUpdatesActive" and will terminate itself)
    this->multipleUpdatesActive = false;
  }


  ////
  //// increment cross section index
  ////
  if (key == GLFW_KEY_TAB) {
    if (this->displayCrossSection) {
      int d = mods == 0 ? 1 : -1;
      int N = Lx + Ly + Lz;
      this->crossSectionIndex = positive_modulo((int) this->crossSectionIndex + d, N);
    }
  }

}

#endif