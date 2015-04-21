#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include <stdio.h>

#include "ion.h"
#include "element.h"

namespace MyTRIM_NS {

struct materialBase {
  double rho;

  double am, az; // average mass and atomic number
  double arho, mu;
  double a, f, epsdg;
  double fd, kd, pmax;

  int tag;
  bool dirty;

  std::vector<elementBase*> element;

  //layerType() { semax = 0.0; sem = 0.0; sez = 0; }
  materialBase( double _rho ) : rho(_rho), tag(-1) { dirty = true; };

  // make sure stoiciometry is normalized, compute averages independent of pka
  void prepare();

  // compute pka dependent averages
  void average( const ionBase *pka );
  double getrstop( const ionBase *pka );

protected:
  double rpstop( int z2, double e );
  double rstop( const ionBase *ion, int z2 );
};

}

#endif
