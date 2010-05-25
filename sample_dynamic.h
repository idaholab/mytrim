#ifndef SAMPLE_DYNAMIC_H
#define SAMPLE_DYNAMIC_H 1

#include "sample_layers.h"
#include <vector>
using namespace std;

struct sampleDynamic : sampleLayers {
  const ionBase *pka;
  vector<bool> layerUpdated;

  virtual void averages( const ionBase *_pka );

  sampleDynamic( float x, float y, float z ): sampleLayers( x, y, z) {};
  virtual materialBase* lookupMaterial( float* pos );
  virtual int lookupLayer( float* pos );
  
  virtual void addAtomsToLayer( int layer, int n, int Z );
};

#endif