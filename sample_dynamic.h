#ifndef SAMPLE_DYNAMIC_H
#define SAMPLE_DYNAMIC_H 1

#include "sample_layers.h"
#include "material.h"
#include <vector>
using namespace std;

struct sampleDynamic : sampleLayers {
  const ionBase *pka;

  virtual void averages( const ionBase *_pka );

  sampleDynamic( double x, double y, double z ): sampleLayers( x, y, z) { bc[0] = CUT; bc[1] = PBC; bc[2] = PBC; };
  virtual materialBase* lookupMaterial( double* pos );

  virtual void addAtomsToLayer( int layer, int n, int Z );
};

#endif