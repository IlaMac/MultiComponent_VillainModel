#ifdef GUI_ENABLED

#include "GUI.h"

using namespace std;


void GUI::onKeyPress (
  const int key,
  const int scancode,
  const int mods
) {

  ////
  //// start/stop updating configuration
  ////
  if (key == GLFW_KEY_ENTER) {
    if ( ! this->multipleUpdatesActive) {
      // multiple update not active
      if (mods == 1) {
      // shift + enter -> multiple updates

        // start to update worm
        this->multipleUpdatesActive = true;
        std::thread( [this] {
          this->tryMultipleUpdates();
        } ).detach();
      } else {
        // enter -> single update
        this->singleUpdate();
      }
    } else {
      // inactivate multiple updates process
      // (the thread listens to "this->multipleUpdatesActive" and will terminate itself)
      this->multipleUpdatesActive = false;
    }
  }


  ////
  //// change component
  ////
  if (key == GLFW_KEY_C) {
    this->actComp = (this->actComp + 1) % (NC + 1);
  }


  ////
  //// display cross section only
  ////
  if (key == GLFW_KEY_X) {
    this->displayCrossSection = ! this->displayCrossSection;
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