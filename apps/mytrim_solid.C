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

int main(int argc, char *argv[])
{
  char fname[200];
  if (argc != 3) // 2
  {
    fprintf(stderr, "syntax:\n%s basename r\n\nCbfactor=1 => 7e-4 bubbles/nm^3\n", argv[0]);
    return 1;
  }

  // seed randomnumber generator from system entropy pool
  FILE *urand = fopen("/dev/random", "r");
  int seed;
  if (fread(&seed, sizeof(int), 1, urand) != 1) return 1;
  fclose(urand);
  r250_init(seed<0 ? -seed : seed); // random generator goes haywire with neg. seed

  // initialize global parameter structure and read data tables from file
  SimconfType * simconf = new SimconfType;
  simconf->fullTraj = false;
  simconf->tmin = 0.2;
  //simconf->tmin = 0.2;

  // initialize sample structure
  SampleSolid *sample = new SampleSolid(200.0, 200.0, 200.0);

  TrimBase *trim = new TrimBase(simconf, sample);

  MaterialBase *material;
  ElementBase *element;

  // UO2
  material = new MaterialBase(simconf, 9.4); // rho
  element = new ElementBase;
  element->_Z = 92; // U
  element->_m = 235.0;
  element->_t = 1.0;
  material->_element.push_back(element);
  element = new ElementBase;
  element->_Z = 16; // O
  element->_m = 32.0;
  element->_t = 2.0;
  material->_element.push_back(element);
  element = new ElementBase;
  element->_Z = 54; // Xe
  element->_m = 131.0;
  element->_t = 0.0024;
  material->_element.push_back(element);
  material->prepare(); // all materials added
  sample->material.push_back(material); // add material to sample

  // create a FIFO for recoils
  std::queue<IonBase*> recoils;

  snprintf(fname, 199, "%s.Erec", argv[1]);
  FILE *erec = fopen(fname, "wt");

  snprintf(fname, 199, "%s.dist", argv[1]);
  FILE *rdist = fopen(fname, "wt");

  IonMDTag *ff1, *pka;

  // 5 fission events
  for (int n = 0; n < 10; n++) // 10 ff
  {
    if (n % 100 == 0) fprintf(stderr, "pka #%d\n", n+1);

    ff1 = new IonMDTag;
    ff1->gen = 0; // generation (0 = PKA)
    ff1->tag = -1;
    ff1->_md = 0;
    ff1->id = simconf->id++;

    ff1->_Z = 53;
    ff1->_m = 127;
    ff1->_E  = 70.0 * 1.0e6;

    Real norm;
    do
    {
      for (int i = 0; i < 3; ++i) ff1->_dir(i) = dr250() - 0.5;
      norm = ff1->_dir.size_sq();
    }
    while (norm <= 0.0001 || norm > 0.25);
    ff1->_dir /= std::sqrt(norm);

    for (int i = 0; i < 3; ++i) ff1->_pos(i) = dr250() * sample->w[i];

    ff1->setEf();
    recoils.push(ff1);

/*
    ff2 = new IonBase(*ff1); // copy constructor
    //ff1->id = simconf->id++;

    // reverse direction
    ff2->_dir = -ff2->_dir;

    ff2->_Z = Z2;
    ff2->_m = A2;
    ff2->_E  = E2 * 1.0e6;

    ff2->setEf();
    recoils.push(ff2);

    fprintf(stderr, "A1=%f Z1=%d (%f MeV)\tA2=%f Z2=%d (%f MeV)\n", A1, Z1, E1, A2, Z2, E2);
*/
    while (!recoils.empty())
    {
      pka = dynamic_cast<IonMDTag*>(recoils.front());
      recoils.pop();
      sample->averages(pka);

      // do ion analysis/processing BEFORE the cascade here

      if (pka->_Z == 54 )
      {
        // mark the first recoil that falls into the MD energy gap with 1 (child generations increase the number)
        if (pka->_E > 200 && pka->_E < 12000 && pka->_md == 0) pka->_md = 1;

        if (pka->gen > 0)
        {
          // output energy and recoil generation
          fprintf(erec, "%f\t%d\t%d\n", pka->_E, pka->gen, pka->_md);
        }

      }

      // follow this ion's trajectory and store recoils

      trim->trim(pka, recoils);

      // do ion analysis/processing AFTER the cascade here

      // pka is Xe
      if (pka->_Z == 54 )
      {
        // output
        //printf("%f %f %f %d\n", pka->_pos(0), pka->_pos(1), pka->_pos(2), pka->tag);
      }

      // done with this recoil
      delete pka;

      // this should rather be done with spawnRecoil returning false
      //if (simconf->primariesOnly) while (!recoils.empty()) { delete recoils.front(); recoils.pop(); };
    }
  }
  fclose(rdist);
  fclose(erec);

  return EXIT_SUCCESS;
}
