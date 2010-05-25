#ifndef SAMPLE_H
#define SAMPLE_H 1

#include <vector>
#include <queue>

#include "ion.h"
#include "material.h"

using namespace std;


struct sampleBase {
  enum sampleBoundary { PBC, INF, CUT }; // periodic, infinitly large, cut off cascades

  vector<materialBase*> material;
  float w[3]; // simulation volume
  sampleBoundary bc[3]; // boundary conditions

  virtual void averages( const ionBase *pka );
  virtual materialBase* lookupMaterial( float* pos ) = 0;
  virtual float rangeMaterial( float* pos, float* dir ) { return 100000.0; };

  sampleBase( float x = 10000.0, float y = 10000.0, float z = 10000.0 );
};

#endif
