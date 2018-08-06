/***************************************************************************
 *   Copyright (C) 2018 by Daniel Schwen   *
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

#include <iostream>
#include <fstream>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>
#include <queue>

#define _USE_MATH_DEFINES
#include <cmath>

#include "simconf.h"
#include "element.h"
#include "material.h"
#include "sample_solid.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"

#include "functions.h"

using namespace MyTRIM_NS;

int
main(int argc, char * argv[])
{
  // initialize global parameter structure and read data tables from file
  SimconfType * simconf = new SimconfType;

  if (argc != 4)
  {
    std::cerr << "syntax:\n"
              << argv[0]
              << " Epka Npka El\nEpka\tPKA energy in eV\nNpka\tNumber of PKA to "
                 "simulate\nEl\tElement symbol for the PKA atom type\n";
    return 1;
  }
  const Real Epka = atof(argv[1]);
  const int Npka = atoi(argv[2]);
  const std::string Elpka(argv[3]);

  // set seed
  int seed;
  char * seedenv = getenv("MYTRIM_SEED");
  if (seedenv)
  {
    // use the number provided in the environment variable MYTRIM_SEED
    seed = atoi(seedenv);
  }
  else
  {
    // seed random number generator from system entropy pool
    // we internally use the libc random function (not r250c, which is not threadsafe)
    FILE * urand = fopen("/dev/random", "r");
    if (fread(&seed, sizeof(int), 1, urand) != 1)
      return 1;
    fclose(urand);
  }
  simconf->seed(seed < 0 ? -seed : seed);

  // initialize sample structure
  auto * sample = new SampleSolid(400.0, 400.0, 400.0);

  // initialize trim engine for the sample
  TrimBase * trim = new TrimBase(simconf, sample);

  MaterialBase * material;
  Element element;

  // site volume
  const Real site_volume = 0.01181;      // nm^3 / Cu_atom
  const Real amu_in_g = 1.660539040e-24; // g / amu
  const Real m = 63.5;                   // amu
  const Real rho = m * amu_in_g / (site_volume * 1e-21);
  std::cerr << "Rho = " << rho << " g/cm^3\n";

  material = new MaterialBase(simconf, 10.0); // rho
  element._Z = 29;                            // Cu
  element._m = m;
  element._t = 1.0;
  material->_element.push_back(element);
  material->prepare();                  // all materials added
  sample->material.push_back(material); // add material to sample

  // create a FIFO for recoils
  std::queue<IonBase *> recoils;

  // distance histogram
  std::vector<int> hist;
  const Real dr = 1.0;

  // find projectile
  Real Mpka = -1, Zpka;
  for (unsigned int i = 0; i < simconf->scoef.size(); ++i)
    if (simconf->scoef[i].sym == Elpka)
    {
      Mpka = simconf->scoef[i].mm1;
      Zpka = i + 1;
      break;
    }
  if (Mpka < 0)
  {
    std::cerr << "Element not found\n";
    return 1;
  }
  std::cerr << "Running " << Elpka << ' ' << Zpka << ' ' << Mpka << '\n';

  Real pos1[3];
  IonMDTag *projectile, *pka;
  for (int n = 0; n < Npka; n++)
  {
    if (n % 10 == 0)
      std::cerr << "event #" << n + 1 << "\n";

    projectile = new IonMDTag;
    projectile->_gen = 0; // generation (0 = PKA)
    projectile->_tag = -1;
    projectile->_md = 0;

    projectile->_Z = Zpka;
    projectile->_m = Mpka;
    projectile->_E = Epka;

    Real norm;
    do
    {
      for (int i = 0; i < 3; ++i)
        projectile->_dir(i) = 2.0 * simconf->drand() - 1.0;
      norm = projectile->_dir.norm_sq();
    } while (norm <= 0.0001 || norm > 1.0);
    projectile->_dir /= std::sqrt(norm);

    // random origin
    for (int i = 0; i < 3; ++i)
      projectile->_pos(i) = simconf->drand() * sample->w[i];

    projectile->setEf();
    recoils.push(projectile);

    while (!recoils.empty())
    {
      pka = dynamic_cast<IonMDTag *>(recoils.front());
      recoils.pop();
      sample->averages(pka);

      for (int i = 0; i < 3; ++i)
        pos1[i] = pka->_pos(i);

      trim->trim(pka, recoils);

      if (pka->_Z == 29)
      {
        Real r2 = 0.0;
        for (int i = 0; i < 3; ++i)
          r2 += (pos1[i] - pka->_pos(i)) * (pos1[i] - pka->_pos(i));
        const std::size_t bin = std::sqrt(r2) / dr;
        if (bin >= hist.size())
          hist.resize(bin + 1);

        hist[bin]++;
      }

      // done with this recoil
      delete pka;
    }
  }

  for (std::size_t i = 0; i < hist.size(); ++i)
    std::cout << i * dr << ' ' << hist[i] / Npka << '\n';

  return EXIT_SUCCESS;
}
