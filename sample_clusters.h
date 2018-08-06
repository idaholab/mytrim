/*
MyTRIM - a three dimensional binary collision Monte Carlo library.
Copyright (C) 2008-2018  Daniel Schwen <daniel@schwen.de>

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

#ifndef MYTRIM_SAMPLE_CLUSTERS_H
#define MYTRIM_SAMPLE_CLUSTERS_H

#include "sample.h"

namespace MyTRIM_NS
{

struct sampleClusters : SampleBase
{

  Real sd, kd[3]; // half the spatial diagonal of a hash block, hash block size
  int *sh, kn[3]; // spatial hash and its dimensions

  int *cl, cn,
      cnm;     // cluster linklist, actual number of clusters (incl. ghosts) and number reserved
  Real * c[4]; // three arrays for cluster x, y, z, r^2 coordinates
  Real cmr;    // maximum cluster radius in the sample

  sampleClusters(Real x = 10000.0, Real y = 10000.0, Real z = 10000.0);

  virtual MaterialBase * lookupMaterial(Point & pos);

  int lookupCluster(Point & pos, Real dr = 0.0);
  void initSpatialhash(int x, int y, int z);
  void addCluster(Real x, Real y, Real z, Real r);
  void addRandomClusters(unsigned int n, Real r, Real dr, SimconfType * simconf);

protected:
  void clearClusters();
  void clearSpatialHash();
  void reallocClusters(int n);
};
}

#endif
