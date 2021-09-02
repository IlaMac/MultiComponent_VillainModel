#ifdef GUI_ENABLED

#ifndef _GUI
#define _GUI

// #include <thread>
// #include <unistd.h>
// #include <chrono>


#include "../../OpenGL/src/Window/Window.h"

#include "../constants.h"
#include "../initialization.h"


using Lattice = std::vector<Node>;


class GUI {
  private:

    Window * window;

    Lattice & lattice;

  public:

  private:

    // the main render loop
    bool renderLoop (const double &);


  public:
    GUI (Lattice &);
    ~GUI ();
};

// the function which creates an instance of the following class and runs it
static inline void runGUI (Lattice & lattice) {
  GUI temp{lattice};
}

#endif

#endif