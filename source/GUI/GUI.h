#ifdef GUI_ENABLED

#ifndef _GUI
#define _GUI

#include <thread>
// #include <unistd.h>
// #include <chrono>
#include <iomanip>
#include <tuple>

#include "../../OpenGL/src/Window/Window.h"

#include "../constants.h"
#include "../initialization.h"
#include "../villain_MC.h"


using Lattice = std::vector<Node>;

template<typename T>
static inline
T positive_modulo(T i, T n) {
  return (n + (i % n)) % n;
}

class GUI {
  private:

    Window * window;

    const Lattice        & lattice;
    struct MC_parameters & MCp;
    struct H_parameters  & Hp;
    double                 beta;
    struct Villain       & villain;


    unsigned actComp = NC;


    ////
    //// controlling multiple updates process
    ////
    bool multipleUpdatesActive = false;
    unsigned long numUpdatesPerFrame = 100; //debug_major ? 1000 : 100000;
    bool blockUpdates = false;
    bool blockRendering = false;


    ////
    //// cross section variables
    ////
    bool displayCrossSection   = false;
    unsigned crossSectionIndex = 0;


    ////
    //// container for velocities etc.
    ////
    std::vector<glm::vec3> vs;
    std::vector<float> maxVs2;

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

    // attempt to perform multiple updates
    void tryMultipleUpdates ();

    // perform a single update/sweep
    void singleUpdate ();

    // the hsv color map
    glm::vec3 colorMapHSV (float) const;

    // compute velocities from phases
    void computeVelocities ();

    // aux
    static
    std::tuple<int, int, int> i2is (const int);

    static
    int is2i (const int, const int, const int);

    static
    int is2i (const int, const int, const int, const int);



  public:
    GUI (const Lattice &,
         struct MC_parameters &,
         struct H_parameters &,
         const double,
         struct Villain &);

    ~GUI ();
};

// the function which creates an instance of the following class and runs it
static inline void runGUI (const Lattice        & lattice,
                           struct MC_parameters & MCp,
                           struct H_parameters  & Hp,
                           const double           beta,
                           struct Villain       & villain
) {
  GUI temp{lattice, MCp, Hp, beta, villain};
}

#endif

#endif