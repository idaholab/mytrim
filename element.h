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

#ifndef MYTRIM_ELEMENT_H
#define MYTRIM_ELEMENT_H

#include "simconf.h"

namespace MyTRIM_NS
{

class Element
{
public:
  Element();

  int _Z;
  Real _m, _t; // mass and relative amount

  Real _Edisp, _Elbind; // displacement energy and lattice binding energy

  // calculated
  Real my, ec, ai, fi;
};
}

#endif
