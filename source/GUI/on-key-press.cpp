#ifdef GUI_ENABLED

#include "GUI.h"

using namespace std;


void GUI::onKeyPress (
  const int key,
  const int scancode,
  const int mods
) {
  // ////
  // //// change component
  // ////
  // if (key == GLFW_KEY_C) {
  //   this->activeComponent = (this->activeComponent + 1) % (numComps + 1);
  // }


  // ////
  // //// draw site placeholders
  // ////
  // if (key == GLFW_KEY_P) {
  //   this->drawSitePlaceholders = ! this->drawSitePlaceholders;
  // }


  // ////
  // //// start/stop updating worm
  // ////
  // if (key == GLFW_KEY_ENTER) {
  //   if ( ! this->multipleUpdatesActive) {
  //     // multiple update not active
  //     if (mods == 1) {
  //     // shift + enter -> multiple updates

  //       // start to update worm
  //       this->multipleUpdatesActive = true;
  //       std::thread( [this] {
  //         this->updateWorm();
  //       } ).detach();
  //     } else {

  //       if ( ! this->sampleOnlyZ) {
  //         // enter -> single update
  //         this->worm.singleUpdate();
  //       } else {
  //         // enter -> update until Z
  //         this->multipleUpdatesActive = true;
  //         std::thread( [this] {
  //           this->updateWorm(true);
  //         } ).detach();
  //       }

  //     }
  //   } else {
  //     // inactivate multiple updates process
  //     // (the thread listens to "this->multipleUpdatesActive" and will terminate itself)
  //     this->multipleUpdatesActive = false;
  //   }
  // }


  // ////
  // //// distinguish fermions
  // ////
  // if (key == GLFW_KEY_F) {
  //   this->distinguishFermions = ! this->distinguishFermions;
  // }


  // ////
  // //// if we should sample only partition functions
  // ////
  // if (key == GLFW_KEY_Z) {
  //   this->sampleOnlyZ = ! this->sampleOnlyZ;
  // }


  // ////
  // //// if we want to outline the head and the tail with a wire frame
  // ////
  // if (key == GLFW_KEY_O) {
  //   this->outlineHeadAndTail = ! this->outlineHeadAndTail;
  // }


  // ////
  // //// toggle simplified drawing
  // ////
  // if (key == GLFW_KEY_L) {
  //   this->simplified = ! this->simplified;
  // }


  // ////
  // //// hide world line
  // ////
  // if (key == GLFW_KEY_H) {
  //   this->hideWormBody = ! this->hideWormBody;
  // }


  // ////
  // //// update parameter from the input parameter file
  // ////
  // if (key == GLFW_KEY_R) {
  //   this->queueUpdateParameters();
  // }


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