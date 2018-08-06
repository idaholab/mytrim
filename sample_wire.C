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

#include <cmath>
#if !defined(__APPLE__)
#include "malloc.h"
#endif

#include "sample_wire.h"

using namespace MyTRIM_NS;

SampleWire::SampleWire(Real x, Real y, Real z) : SampleBase(x, y, z)
{
  bc[0] = CUT;
  bc[1] = CUT;
}

// look if we are within dr of the wire axis
MaterialBase *
SampleWire::lookupMaterial(Point & pos)
{
  Real x = (pos(0) / w[0]) * 2.0 - 1.0;
  Real y = (pos(1) / w[1]) * 2.0 - 1.0;

  if ((x * x + y * y) > 1.0)
    return 0;
  return material[0];
}
