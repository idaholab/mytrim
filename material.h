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

struct materialBase {
  double rho;

  // set in prepare
  double am, az; // average mass and atomic number
  double arho;

  // set in average
  double mu;
  double a, f, epsdg;
  double fd, kd;

  double pmax;

  int tag;
  bool dirty;

  std::vector<elementBase*> element;

  materialBase(simconfType * simconf_, double rho_);

  // make sure stoiciometry is normalized, compute averages independent of pka
  void prepare();

  // compute pka dependent averages
  void average(const ionBase *pka);
  double getrstop(const ionBase *pka);

  virtual elementBase * getElement(unsigned int nn) { return element[nn]; }

protected:
  double rpstop(int z2, double e);
  double rstop(const ionBase *ion, int z2);

  simconfType * simconf;
};

}

#endif
