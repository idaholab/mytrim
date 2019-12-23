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

#ifndef MYTRIM_SAMPLE_H
#define MYTRIM_SAMPLE_H

#include <vector>
#include <queue>

#include "ion.h"
#include "material.h"

namespace MyTRIM_NS
{

class SampleBase
{
public:
  SampleBase(Real x = 10000.0, Real y = 10000.0, Real z = 10000.0);

  virtual void averages(const IonBase * pka);
  virtual MaterialBase * lookupMaterial(Point & pos) = 0;
  virtual Real rangeMaterial(Point & pos, Point & dir);

  /// materials in the sample
  std::vector<MaterialBase *> material;

  /// simulation volume
  Real w[3];

  /// periodic, infinitly large, cut off cascades
  enum sampleBoundary
  {
    PBC,
    INF,
    CUT
  };

  /// boundary conditions
  sampleBoundary bc[3];
};
} // namespace MyTRIM_NS

#endif
