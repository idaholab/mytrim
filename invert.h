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

#ifndef MYTRIM_INVERT_H
#define MYTRIM_INVERT_H

#include "simconf.h"
#include <cmath>

namespace MyTRIM_NS {

class inverter
{
protected :
  Real maxx, maxf, tol;

public :
  virtual Real f(Real x) = 0;
  Real x(Real f);

  inverter() { maxx = 0.0; }
};


class massInverter : public inverter
{
public:
  virtual Real f(Real x);

  massInverter() { maxx = 235.0; tol = 1e-7; maxf = f(maxx); }
};

class energyInverter : public inverter
{
 Real A;
public:
  virtual Real f(Real x);

  energyInverter() { maxx = 186.98; tol = 1e-7; setMass(100.0); }
  void setMass(Real _A) {
    A = _A;
    maxf = f(maxx);
  }
};

}

#endif
