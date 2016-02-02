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

#ifndef MYTRIM_SIMCONF_H
#define MYTRIM_SIMCONF_H

#ifdef MYTRIM_ENABLED
// building from within MOOSE/Magpie
#include "MooseError.h"
#include "MooseTypes.h"
#include "libmesh/point.h"
namespace MyTRIM_NS {
  const Real drm = Real(RAND_MAX)+1.0;
  inline Real dr250() { return Real(rand())/drm; }
  inline void r250_init(int s) { srand(s); }
}
#else
// building standalone (for Travis CI tests)
typedef double Real;
#include "point.h"
#include "r250.h"
#endif

#include <stdlib.h>

namespace MyTRIM_NS {

// ZBL coefficients a,b,c,d for all element pairs from Z=1..92
struct scoefLine {
  char sym[3], name[30];
  Real mm1, m1, mnat, rho, atrho, vfermi, heat, lfctr;
};

struct simconfType {
  Real ed, alfa, alpha, tmin, tau, da, cw;
  int id;

  // tables from files
  scoefLine scoef[92];
  Real pcoef[92][8];
  Real snuc[92][92][4];

  bool fullTraj;

  // tally electronic and nuclear stopping losses
  Real EelTotal;
  Real EnucTotal;

  // statistics of the simulation run
  int vacancies_created;

  simconfType(Real _alfa = 0.0);
private:
  void read_scoef();
  void read_snuc();

  void fileReadError(const char *);
};

}

#endif
