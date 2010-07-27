#ifndef ELEMENT_H
#define ELEMENT_H 1

struct elementBase {
  int z;
  float m, t; // mass and relative amount

  float Edisp, Elbind; // displacement energy and lattice binding energy

  // calculated
  float my, ec, ai, fi;

  elementBase();
};


#endif