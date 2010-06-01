#ifndef SAMPLE_LAYERS_H
#define SAMPLE_LAYERS_H 1

#include "sample.h"
#include <vector>
using namespace std;

struct sampleLayers : sampleBase {
  vector<double> layerThickness;

  sampleLayers( float x, float y, float z ): sampleBase( x, y, z) {};
  virtual int lookupLayer( float* pos );
  virtual materialBase* lookupMaterial( float* pos );
  virtual float rangeMaterial( float* pos, float* dir );
};

#endif