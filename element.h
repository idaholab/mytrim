#ifndef MYTRIM_ELEMENT_H
#define MYTRIM_ELEMENT_H

namespace MyTRIM_NS {

struct elementBase {
  int z;
  double m, t; // mass and relative amount

  double Edisp, Elbind; // displacement energy and lattice binding energy

  // calculated
  double my, ec, ai, fi;

  elementBase();
};

}

#endif
