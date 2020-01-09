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

#ifndef MYTRIM_INVERT_H
#define MYTRIM_INVERT_H

#include "simconf.h"
#include <cmath>

namespace MyTRIM_NS
{

class Inverter
{
public:
  Inverter() : maxx(0.0), maxf(0.0), tol(1e-13) {}

  /// Evaluate inverse of function (iteratively)
  Real x(Real f) const;

protected:
  /// Evaluate function
  virtual Real f(Real x) const = 0;

  Real maxx, maxf, tol;
};

class MassInverter : public Inverter
{
public:
  MassInverter()
  {
    maxx = 235.0;
    tol = 1e-7;
    maxf = f(maxx);
  }

protected:
  virtual Real f(Real x) const;
};

class EnergyInverter : public Inverter
{
public:
  EnergyInverter()
  {
    maxx = 186.98;
    tol = 1e-7;
    setMass(100.0);
  }

  void setMass(Real A)
  {
    _A = A;
    maxf = f(maxx);
  }

protected:
  virtual Real f(Real x) const;

private:
  Real _A;
};
} // namespace MyTRIM_NS

#endif
