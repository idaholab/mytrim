#ifndef SAMPLE_LAYERS_H
#define SAMPLE_LAYERS_H

#include "sample.h"
#include <vector>

namespace MyTRIM_NS {

struct sampleLayers : sampleBase {
  std::vector<double> layerThickness;

  sampleLayers( double x, double y, double z ): sampleBase( x, y, z) {};
  virtual int lookupLayer( double* pos );
  virtual materialBase* lookupMaterial( double* pos );
  virtual double rangeMaterial( double* pos, double* dir );
};

}

#endif
