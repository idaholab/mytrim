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

#ifndef MYTRIM_MATERIAL_H
#define MYTRIM_MATERIAL_H

#include <vector>
#include <stdio.h>

#include "ion.h"
#include "element.h"
#include "simconf.h"

namespace MyTRIM_NS {

class MaterialBase
{
public:
  MaterialBase(SimconfType * simconf, Real rho);

  /// make sure stoiciometry is normalized, compute averages independent of pka
  void prepare();

  /// compute pka dependent averages
  void average(const IonBase *pka);
  Real getrstop(const IonBase *pka);

  virtual ElementBase * getElement(unsigned int nn) { return _element[nn]; }

  Real _rho;

  // set in prepare
  Real _am, _az; // average mass and atomic number
  Real _arho;

  // set in average
  Real mu;
  Real a, f, epsdg;
  Real fd, kd;

  Real pmax;

  int tag;
  bool _dirty;

  std::vector<ElementBase*> _element;

protected:
  Real rpstop(int z2, Real e);
  Real rstop(const IonBase *ion, int z2);

  SimconfType * _simconf;
};

}

#endif
