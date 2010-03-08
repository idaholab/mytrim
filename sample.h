#ifndef SAMPLE_H
#define SAMPLE_H 1

#include <vector>
#include <queue>

#include "ion.h"
#include "material.h"

using namespace std;

enum sampleBoundary { PBC, INF, CUT }; // periodic, infinitly large, cut off cascades

struct sampleBase {
  vector<materialBase*> material;
  float w[3]; // simulation volume
  sampleBoundary bc[3]; // boundary conditions

  void averages( const ionBase *pka );
  virtual materialBase* lookupMaterial( float* pos ) = 0;

  sampleBase( float x = 10000.0, float y = 10000.0, float z = 10000.0 );
};

#endif
