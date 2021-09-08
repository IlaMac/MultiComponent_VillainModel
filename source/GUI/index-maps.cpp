#ifdef GUI_ENABLED

#include "GUI.h"

using namespace std;

// inline
tuple<int, int, int>
GUI::i2is (
  const int i
) {
  cout << "implement i2is()" << endl;
  exit(1);
  return {Lx, Ly, Lz};
}


// inline
int
GUI::is2i (
  const int ix,
  const int iy,
  const int iz
) {
  return ix + Lx * (iy + iz * Ly);
}


int
GUI::is2i (
  const int a,
  const int ix,
  const int iy,
  const int iz
) {
  return a*Lx*Ly*Lz + is2i(ix, iy, iz);
}

#endif