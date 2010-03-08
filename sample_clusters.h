#ifndef SAMPLE_CLUSTERS_H
#define SAMPLE_CLUSTERS_H 1

#include "sample.h"

struct sampleClusters : sampleBase {

  float sd, kd[3]; // half the spatial diagonal of a hash block, hash block size
  int *sh, kn[3]; // spatial hash and its dimensions

  int *cl, cn, cnm; // cluster linklist, actual number of clusters (incl. ghosts) and number reserved
  float *c[4]; // three arrays for cluster x,y,z,r^2 coordinates
  float cmr; // maximum cluster radius in the sample

  sampleClusters( float x = 10000.0, float y = 10000.0, float z = 10000.0 );

  virtual materialBase* lookupMaterial( float* pos );

  int lookupCluster( float* pos, float dr = 0.0 );
  void initSpatialhash( int x, int y, int z );
  void clearSpatialHash();
  void clearClusters();
  void addCluster( float x, float y, float z, float r );
  void addRandomClusters( int n, float r, float dr = 0.0 );
protected:
  void reallocClusters( int n );
};

#endif