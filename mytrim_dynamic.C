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
#include "sample_dynamic.h"
#include "ion.h"
#include "trim.h"

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
  fread(&seed, sizeof(int), 1, urand);
  fclose(urand);
  seed = 1;
  r250_init(seed<0 ? -seed : seed); // random generator goes haywire with neg. seed

  // initialize global parameter structure and read data tables from file
  simconfType * simconf = new simconfType;
  simconf->fullTraj = false;
  simconf->tmin = 0.2;
  //simconf->tmin = 0.2;

  // initialize sample structure
  double sx, sy, sz;
  std::cin >> sx >> sy >> sz;
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  std::cout << "SS " << sx << ' ' << sy << ' ' << sz << std::endl;

  sampleDynamic *sample = new sampleDynamic(simconf, sx, sy, sz);
  trimBase *trim = new trimBase(simconf, sample);

  // Read Materials description from stdin
  int nlayer;
  double lthick, lrho, nelem;
  std::string lename;
  std::cin >> nlayer;
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  std::cout << "n_layers=" << nlayer << std::endl;

  materialBase *material;
  elementBase *element;

  for (int i = 0; i < nlayer; i++)
  {
    std::cin >> lename >> lthick >> lrho >> nelem;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout << "Layer: " << lename << "  d=" << lthick << "Ang  rho="
         << lrho << "g/ccm  n_elements=" << nelem << std::endl;

    material = new materialBase(simconf, lrho); // rho

    for (int j = 0; j < nelem; j++)
    {
      element = new elementBase;
      std::cin >> lename >> element->z >> element->m >> element->t;
      std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      std::cout << "  Element: " << lename << "  Z=" << element->z
           << "  m=" << element->m << "  fraction=" << element->t << std::endl;
      material->element.push_back(element);
    }
    std::cout << material << std::endl;

    material->prepare(); // all elements added
    sample->material.push_back(material); // add material to sample
    sample->layerThickness.push_back(lthick);
  }

  // create a FIFO for recoils
  std::queue<ionBase*> recoils;

  double norm;
  double jmp = 2.7; // diffusion jump distance
  int jumps;
  double dif[3];

  // double A = 74.0, E = 1.0e5; int Z = 36; // 100keV Kr
  // double A = 131.0, E = 3.0e4; int Z = 54; // 30keV Xe
  double A = 131.0, E = 1.0e4; int Z = 54; // 30keV Xe

  snprintf(fname, 199, "%s.Erec", argv[1]);
  FILE *erec = fopen(fname, "wt");

  snprintf(fname, 199, "%s.dist", argv[1]);
  FILE *rdist = fopen(fname, "wt");

  double pos1[3];

  ionBase *ff1, *ff2, *pka;
  int id = 1;
  int layer1, layer2;

  int nrec = 0;
  double sum_r2, opos[3];

  // 1000 PKA
  for (int n = 0; n < 35000; n++)
  {
    if (n % 100 == 0) fprintf(stderr, "pka #%d\n", n+1);

    ff1 = new ionBase;
    ff1->gen = 0; // generation (0 = PKA)
    ff1->tag = -1;
    ff1->id = simconf->id++;

    ff1->z1 = Z;
    ff1->m1 = A;
    ff1->e  = E * (3-((n*3)/35000));

    ff1->dir[0] = 1;
    ff1->dir[1] = 0;
    ff1->dir[2] = 0;

    ff1->pos[0] = 0;
    ff1->pos[1] = sample->w[1] / 2.0;
    ff1->pos[2] = sample->w[2] / 2.0;

    ff1->set_ef();
    recoils.push(ff1);

    while (!recoils.empty())
    {
      pka = recoils.front();
      recoils.pop();

      sample->averages(pka);

      //
      // do ion analysis/processing BEFORE the cascade here
      //

      // get layer of origin
      layer1 = sample->lookupLayer(pka->pos);

      // remove from source layer
      if (pka->gen > 0)
        sample->addAtomsToLayer(layer1, -1, pka->z1);

      //fprintf(erec, "%f\t%d\t%d\n", pka->e, pka->gen, pka->z1);
      //for (int i = 0; i < 3; i++)
      // opos[i] = pka->pos[i];

      //
      // follow this ion's trajectory and store recoils
      //
      trim->trim(pka, recoils);

      // add to destination layer
      layer2 = sample->lookupLayer(pka->pos);
      if (pka->pos[0] > 0)
        sample->addAtomsToLayer(layer2, 1, pka->z1);

      //
      // do ion analysis/processing AFTER the cascade here
      //

      // done with this recoil
      delete pka;
    }
  }
  fclose(rdist);
  fclose(erec);

  for (int i = 0; i < sample->material.size(); i++)
  {
    std::cout << sample->layerThickness[i] << ' ';
    for (int j = 0; j < sample->material[i]->element.size(); j++)
    {
      std::cout << sample->material[i]->element[j]->z << ' ' << sample->material[i]->element[j]->t << ' ';
    }
    std::cout << std::endl;
  }

  std::cout << "n=" << nrec << " sum_r2=" << sum_r2 << std::endl;
  return EXIT_SUCCESS;
}
