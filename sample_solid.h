#ifndef SAMPLE_SOLID_H
#define SAMPLE_SOLID_H 1

#include "sample.h"

struct sampleSolid : sampleBase {
  sampleSolid( double x, double y, double z ): sampleBase( x, y, z) {};
  virtual materialBase* lookupMaterial( double* pos );
};

#endif