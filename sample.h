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
  double w[3]; // simulation volume
  sampleBoundary bc[3]; // boundary conditions

  virtual void averages( const ionBase *pka );
  virtual materialBase* lookupMaterial( double* pos ) = 0;
  virtual double rangeMaterial( double* pos, double* dir ) { return 100000.0; };

  sampleBase( double x = 10000.0, double y = 10000.0, double z = 10000.0 );
};

#endif
