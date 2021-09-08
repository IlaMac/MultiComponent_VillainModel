#ifdef GUI_ENABLED

#include "GUI.h"

using namespace std;


void GUI::tryMultipleUpdates () {

  unsigned long prevFrameNum = this->window->getFrameNumber();

  ////
  //// loop until it is told to shut down
  ////
  while (this->multipleUpdatesActive) {

    ////
    //// wait for new frame and ensure we are allowed to update
    ////
    if (prevFrameNum != this->window->getFrameNumber() && ! this->blockUpdates) {
      prevFrameNum = this->window->getFrameNumber();

      // disable rendering
      this->blockRendering = true;

      // update normally
      for (unsigned long i = 0; i < this->numUpdatesPerFrame; i++) {
        this->singleUpdate();
      }

      // enabled rendering
      this->blockRendering = false;
    }

  }
}

#endif