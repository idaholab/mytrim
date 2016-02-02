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
#include <iostream>
#include <string>
#include <limits>

#include "simconf.h"
#include "element.h"
#include "material.h"
#include "sample_layers.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"

#include "functions.h"

using namespace MyTRIM_NS;

int main(int argc, char *argv[])
{
  char fname[200];
  if (argc != 2) // 2
  {
    std::cerr << "syntax:\n" << argv[0] << " basename" << std::endl;
    return 1;
  }

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
  Real sx, sy, sz;
  std::cin >> sx >> sy >> sz;
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  std::cout << "SS " << sx << ' ' << sy << ' ' << sz << std::endl;

  const int nmax = 10000;
  std::cout << "NN " << nmax << " PKAs" << std::endl;

  SampleLayers *sample = new SampleLayers(sx, sy, sz);
  //trimBase *trim = new trimBase(sample);
  trimBase *trim = new trimRecoils(simconf, sample);

  // Read Materials description from stdin
  int nlayer;
  Real lthick, lrho, nelem;
  std::string lename;
  std::cin >> nlayer;
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  std::cout << "n_layers=" << nlayer << std::endl;

  MaterialBase *material;
  ElementBase *element;
  for (int i = 0; i < nlayer; i++)
  {
    std::cin >> lename >> lthick >> lrho >> nelem;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout << "Layer: " << lename << "  d=" << lthick << "Ang  rho="
         << lrho << "g/ccm  n_elements=" << nelem << std::endl;

    material = new MaterialBase(simconf, lrho); // rho

    for (int j = 0; j < nelem; j++)
    {
      element = new ElementBase;
      std::cin >> lename >> element->_Z >> element->_m >> element->_t;
      std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      std::cout << "  Element: " << lename << "  Z=" << element->_Z
           << "  m=" << element->_m << "  fraction=" << element->_t << std::endl;
      material->element.push_back(element);
    }

    material->prepare(); // all elements added
    sample->material.push_back(material); // add material to sample
    sample->layerThickness.push_back(lthick);
  }

  // create a FIFO for recoils
  std::queue<IonBase*> recoils;

  //Real A = 74.0, E = 1.0e5; int Z = 36; // 100keV Kr
  Real A = 131.0, E = 5.0e5; int Z = 54; // 500keV Xe

  snprintf(fname, 199, "%s.Erec", argv[1]);
  FILE *erec = fopen(fname, "wt");

  snprintf(fname, 199, "%s.dist", argv[1]);
  FILE *rdist = fopen(fname, "wt");

  IonBase *ff1, *pka;
  int nrec = 0;
  Real sum_r2 = 0.0;
  Point opos;

  // 1000 PKA
  for (int n = 0; n < nmax; n++)
  {
    if (n % 100 == 0) fprintf(stderr, "pka #%d\n", n+1);

    ff1 = new IonBase;
    ff1->gen = 0; // generation (0 = PKA)
    ff1->tag = -1;
    ff1->id = simconf->id++;

    ff1->_Z = Z;
    ff1->_m = A;
    ff1->e  = E;

    ff1->dir(0) = 1;
    ff1->dir(1) = 0;
    ff1->dir(2) = 0;

    ff1->pos(0) = 0;
    ff1->pos(1) = sample->w[1] / 2.0;
    ff1->pos(2) = sample->w[2] / 2.0;

    ff1->setEf();
    recoils.push(ff1);

    while (!recoils.empty())
    {
      pka = recoils.front();
      recoils.pop();
      sample->averages(pka);

      // do ion analysis/processing BEFORE the cascade here

      //fprintf(erec, "%f\t%d\t%d\n", pka->e, pka->gen, pka->_Z);

      opos = pka->pos;

      // follow this ion's trajectory and store recoils
      //if (pka->_Z == 29 || pka->_Z == Z)
      //if (pka->_Z == 29 || pka->_Z == Z)
      trim->trim(pka, recoils);

      // do ion analysis/processing AFTER the cascade here
      if (pka->_Z != Z)
      {
        sum_r2 += (opos - pka->pos).size_sq();
        nrec++;
      }

      // pka is O or Ag
      //if (pka->_Z == 29 && pka->pos(0) >= 500.0)
      if (pka->_Z == 29)
      {
        // output
        printf("RP %f %d %d\n", pka->pos(0), n,  pka->gen);
      }

      // done with this recoil
      delete pka;

      // this should rather be done with spawnRecoil returning false
      //if (simconf->primariesOnly) while (!recoils.empty()) { delete recoils.front(); recoils.pop(); };
    }
  }
  fclose(rdist);
  fclose(erec);

  std::cout << "n=" << nrec << " sum_r2=" << sum_r2 << std::endl;
  return EXIT_SUCCESS;
}
