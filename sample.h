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

#ifndef MYTRIM_SAMPLE_H
#define MYTRIM_SAMPLE_H

#include <vector>
#include <queue>

#include "ion.h"
#include "material.h"

namespace MyTRIM_NS {

struct sampleBase {
  enum sampleBoundary { PBC, INF, CUT }; // periodic, infinitly large, cut off cascades

  std::vector<materialBase*> material;
  Real w[3]; // simulation volume
  sampleBoundary bc[3]; // boundary conditions

  virtual void averages(const ionBase  * pka);
  virtual materialBase* lookupMaterial(Real * pos) = 0;
  virtual Real rangeMaterial(Real * /* pos */, Real * /* dir */) { return 100000.0; };

  sampleBase(Real x = 10000.0, Real y = 10000.0, Real z = 10000.0);
};

}

#endif
