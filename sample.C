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

#include "sample.h"

using namespace MyTRIM_NS;

SampleBase::SampleBase(Real x, Real y, Real z)
{
  w[0] = x;
  w[1] = y;
  w[2] = z;
  for (unsigned int i = 0; i < 3; ++i)
    bc[i] = PBC;
}

void
SampleBase::averages(const IonBase * pka)
{
  for (unsigned int i = 0; i < material.size(); ++i)
    material[i]->average(pka);
}

Real
SampleBase::rangeMaterial(Point & /* pos */, Point & /* dir */)
{
  return 100000.0;
}
