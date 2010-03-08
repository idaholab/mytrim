#ifndef SAMPLE_SOLID_H
#define SAMPLE_SOLID_H 1

#include "sample.h"

struct sampleSolid : sampleBase {
  sampleSolid( float x, float y, float z ): sampleBase( x, y, z) {};
  virtual materialBase* lookupMaterial( float* pos );
};

#endif