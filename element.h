#ifndef ELEMENT_H
#define ELEMENT_H 1

struct elementBase {
  int z;
  double m, t; // mass and relative amount

  double Edisp, Elbind; // displacement energy and lattice binding energy

  // calculated
  double my, ec, ai, fi;

  elementBase();
};


#endif