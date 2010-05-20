#ifndef SAMPLE_SOLID_H
#define SAMPLE_SOLID_H 1

#include "sample.h"
#include <vector>
using namespace std;

struct sampleLayers : sampleBase {
  vector<double> layerThickness;

  sampleLayers( float x, float y, float z ): sampleBase( x, y, z) {};
  virtual materialBase* lookupMaterial( float* pos );
};

#endif