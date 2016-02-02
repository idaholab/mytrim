#ifndef MYTRIM_ELEMENT_H
#define MYTRIM_ELEMENT_H

#include "simconf.h"

namespace MyTRIM_NS {

class ElementBase
{
public:
  ElementBase();

  int _Z;
  Real _m, _t; // mass and relative amount

  Real _Edisp, _Elbind; // displacement energy and lattice binding energy

  // calculated
  Real my, ec, ai, fi;
};

}

#endif
