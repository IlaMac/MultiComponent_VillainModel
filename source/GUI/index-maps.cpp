#ifdef GUI_ENABLED

#include "GUI.h"

using namespace std;

tuple<int, int, int>
GUI::i2is (
  const int i
) {
  throwError() << "implement me!" << endl;

  return {Lx, Ly, Lz};
}


int
GUI::is2i (
  const int ix,
  const int iy,
  const int iz
) {
  // return ix + Lx * (iy + iz * Ly);
  return iz*Lx*Ly + iy*Lx + ix;
}


int
GUI::is2i (
  const int a,
  const int ix,
  const int iy,
  const int iz
) {
  // return a*Lx*Ly*Lz + is2i(ix, iy, iz);
  // return ix + Lx * (iy + iz * Ly);
  return iz*Lx*Ly*NC + iy*Lx*NC + ix*NC + a;
}


tuple<unsigned, unsigned>
GUI::i2iaib (
  const unsigned I
) {
  // linear index -> upper triangular (ia, ib) index
  const unsigned ia = NC - 2 - floor(sqrt(-8*I + 4*NC*(NC-1)-7)/2.0 - 0.5);
  const unsigned ib = I + ia + 1 - NC*(NC-1)/2 + (NC-ia)*((NC-ia)-1)/2;
  return {ia, ib};
}

unsigned
GUI::iaib2i (
  const unsigned ia,
  const unsigned ib
) {
  // upper triangular (ia, ib) index -> linear index
  return (NC*(NC-1)/2) - (NC-ia)*((NC-ia)-1)/2 + ib - ia - 1;
}




#endif