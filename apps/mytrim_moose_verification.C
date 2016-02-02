/***************************************************************************
 *   Copyright (C) 2008 by Daniel Schwen   *
 *   daniel@schwen.de   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <queue>

#include "simconf.h"
#include "element.h"
#include "material.h"
#include "sample_solid.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"

#include "functions.h"

using namespace MyTRIM_NS;

int main(int, char **)
{
  // seed randomnumber generator from system entropy pool
  FILE *urand = fopen("/dev/random", "r");
  int seed;
  if (fread(&seed, sizeof(int), 1, urand) != 1) return 1;
  fclose(urand);
  r250_init(seed<0 ? -seed : seed); // random generator goes haywire with neg. seed

  // initialize global parameter structure and read data tables from file
  simconfType * simconf = new simconfType;
  simconf->fullTraj = false;
  simconf->tmin = 0.2;
  //simconf->tmin = 0.2;

  // initialize sample structure
  sampleSolid *sample = new sampleSolid(200.0, 200.0, 200.0);

  trimBase *trim = new trimBase(simconf, sample);

  sample->bc[0] = sampleBase::CUT; // no PBC in x (just clusterless matrix)
  sample->bc[1] = sampleBase::CUT; // no PBC in x (just clusterless matrix)
  sample->bc[2] = sampleBase::CUT; // no PBC in x (just clusterless matrix)

  MaterialBase *material;
  ElementBase *element;

  material = new MaterialBase(simconf, 1.0); // rho
  element = new ElementBase;
  element->_Z = 20;
  element->_m = 40.0;
  element->_t = 1.0;
  material->element.push_back(element);
  material->prepare(); // all materials added
  sample->material.push_back(material); // add material to sample

  // create a FIFO for recoils
  std::queue<IonBase*> recoils;

  // create a bunch of ions
  MyTRIM_NS::IonBase * pka;
  for (unsigned int i = 0; i < 1000; ++i)
  {
    pka = new MyTRIM_NS::IonBase;
    pka->gen = 0;  // generation (0 = PKA)
    pka->tag = 0; // tag holds the element type
    pka->_Z = 20;
    pka->_m = 40;
    pka->e  = 300;

    pka->dir(0) = 0.0;
    pka->dir(1) = 1.0;
    pka->dir(2) = 0.0;

    pka->pos(0) = 100.0;
    pka->pos(1) = 0.01;
    pka->pos(2) = 100.0;

    pka->setEf();
    recoils.push(pka);
  }

  while (!recoils.empty())
  {
    pka = recoils.front();
    recoils.pop();
    sample->averages(pka);

    pka->pos(2) = 100.0;
    pka->dir(2) = 0.0;

    trim->trim(pka, recoils);

    printf("%f %f %f\n", pka->pos(0), pka->pos(1), pka->pos(2));

    // done with this recoil
    delete pka;
  }

  return EXIT_SUCCESS;
}
