#ifndef SAMPLE_WIRE_H
#define SAMPLE_WIRE_H 1

#include "sample.h"

struct sampleWire : sampleBase {
  sampleWire( float x = 10000.0, float y = 10000.0, float z = 10000.0 );

  virtual materialBase* lookupMaterial( float* pos );
};

#endif