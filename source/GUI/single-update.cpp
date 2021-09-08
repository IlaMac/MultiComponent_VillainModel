#ifdef GUI_ENABLED

#include "GUI.h"

using namespace std;


void GUI::singleUpdate () {
  metropolis_villain(this->lattice, this->MCp, this->Hp, this->beta, this->villain);
}

#endif