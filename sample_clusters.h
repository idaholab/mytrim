/*
MyTRIM - a three dimensional binary collision Monte Carlo library.
Copyright (C) 2008-2015  Daniel Schwen <daniel@schwen.de>

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA
*/

#ifndef SAMPLE_CLUSTERS_H
#define SAMPLE_CLUSTERS_H

#include "sample.h"

namespace MyTRIM_NS {

struct sampleClusters : sampleBase {

  double sd, kd[3]; // half the spatial diagonal of a hash block, hash block size
  int *sh, kn[3]; // spatial hash and its dimensions

  int *cl, cn, cnm; // cluster linklist, actual number of clusters (incl. ghosts) and number reserved
  double *c[4]; // three arrays for cluster x,y,z,r^2 coordinates
  double cmr; // maximum cluster radius in the sample

  sampleClusters( double x = 10000.0, double y = 10000.0, double z = 10000.0 );

  virtual materialBase* lookupMaterial( double* pos );

  int lookupCluster( double* pos, double dr = 0.0 );
  void initSpatialhash( int x, int y, int z );
  void clearSpatialHash();
  void clearClusters();
  void addCluster( double x, double y, double z, double r );
  void addRandomClusters( int n, double r, double dr = 0.0 );
protected:
  void reallocClusters( int n );
};

}

#endif
