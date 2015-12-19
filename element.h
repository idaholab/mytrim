#ifndef MYTRIM_ELEMENT_H
#define MYTRIM_ELEMENT_H

#include "simconf.h"

namespace MyTRIM_NS {

struct elementBase {
  int z;
  Real m, t; // mass and relative amount

  Real Edisp, Elbind; // displacement energy and lattice binding energy

  // calculated
  Real my, ec, ai, fi;

  elementBase();
};

}

#endif
