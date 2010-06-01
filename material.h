#ifndef MATERIAL_H
#define MATERIAL_H 1

#include <vector>
#include <stdio.h>

#include "ion.h"
#include "element.h"

using namespace std;

struct materialBase {
  float rho;

  float am, az; // average mass and atomic number
  float arho, mu;
  float a, f, epsdg;
  float fd, kd, pmax;

  int tag;
  bool dirty;

  vector<elementBase*> element;

  //layerType() { semax = 0.0; sem = 0.0; sez = 0; }
  materialBase( float _rho ) : rho(_rho) { dirty = true; };

  // make sure stoiciometry is normalized, compute averages independent of pka
  void prepare();

  // compute pka dependent averages 
  void average( const ionBase *pka );
  float getrstop( const ionBase *pka );

protected:
  float rpstop( int z2, float e );
  float rstop( const ionBase *ion, int z2 );
};

#endif
