#ifdef GUI_ENABLED

#ifndef _GUI
#define _GUI

// #include <thread>
// #include <unistd.h>
// #include <chrono>
#include <iomanip>

#include "../../OpenGL/src/Window/Window.h"

#include "../constants.h"
#include "../initialization.h"


using Lattice = std::vector<Node>;

template<typename T>
static inline
T positive_modulo(T i, T n) {
  return (n + (i % n)) % n;
}

class GUI {
  private:

    Window * window;

    Lattice & lattice;

    ////
    //// cross section variables
    ////
    bool displayCrossSection   = false;
    unsigned crossSectionIndex = 0;

  public:

  private:

    // the main render loop
    bool renderLoop (const double &);

    // event listeners
    void onKeyPress (const int,
                     const int,
                     const int);
    void onKeyRepeat (const int,
                      const int,
                      const int);


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