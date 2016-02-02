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
#include "sample_clusters.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"

#include "functions.h"

using namespace MyTRIM_NS;

int main(int argc, char *argv[])
{
  char fname[200];
  if (argc != 5) {
    std::cerr << "syntax:\n"
              << argv[0] << " basename r Cbfactor Nev\n\n"
              << "r Bubble radius in Ang\n"
              << "Cbfactor=1 => 7e-4 bubbles/nm^3\n"
              << "Nev  number of fission events (two fragemnts each)\n";
    return 1;
  }

  // run mode
  enum RunMode { PLAIN, PHONONS, DEFECTS };
  RunMode mode = PLAIN;
  //RunMode mode = PHONONS;

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
    FILE *urand = fopen("/dev/random", "r");
    if (fread(&seed, sizeof(int), 1, urand) != 1) return 1;
    fclose(urand);
  }
  r250_init(seed<0 ? -seed : seed);

  // initialize global parameter structure and read data tables from file
  simconfType * simconf = new simconfType;

  // initialize sample structure
  sampleClusters *sample = new sampleClusters(400.0, 400.0, 400.0);

  // initialize trim engine for the sample
  trimBase *trim;
  std::ofstream auxout;
  std::stringstream auxoutname;
  switch (mode) {
    case PLAIN:
      trim = new trimBase(simconf, sample);
      break;

    case PHONONS:
      auxoutname << argv[1] << ".phonons";
      auxout.open(auxoutname.str().c_str());
      trim = new trimPhononOut(simconf, sample, auxout);
      break;

    case DEFECTS:
      auxoutname << argv[1] << ".defects";
      auxout.open(auxoutname.str().c_str());
      trim = new trimDefectLog(simconf, sample, auxout);
      break;

    default:
      return 1;
  }

  //Real r = 10.0;
  Real r = atof(argv[2]); //10.0;
  Real Cbf = atof(argv[3]);
  int Nev = atoi(argv[4]);

  //sample->bc[0] = CUT; // no pbc in x dir
  sample->initSpatialhash(int(sample->w[0] / r) - 1,
                           int(sample->w[1] / r) - 1,
                           int(sample->w[2] / r) - 1);


  // Real atp = 0.1; // 10at% Mo 90at%Cu
  Real v_sam = sample->w[0] * sample->w[1] * sample->w[2];
  //Real v_cl = 4.0/3.0 * M_PI * cub(r);
  int n_cl; // = atp * scoef[29-1].atrho * v_sam / (v_cl * ((1.0 - atp) * scoef[42-1].atrho + atp * scoef[29-1].atrho));

  n_cl = v_sam * 7.0e-7 * Cbf ; // Ola06 7e-4/nm^3
  std::cerr << "adding " << n_cl << " clusters...\n";

  // cluster surfaces must be at least 25.0 Ang apart
  sample->addRandomClusters(n_cl, r, 25.0);

  // write cluster coords with tag numbers
  snprintf(fname, 199, "%s.clcoor", argv[1]);
  FILE *ccf = fopen(fname, "wt");
  for (int i = 0; i < sample->cn; i++)
    fprintf(ccf, "%f %f %f %f %d\n", sample->c[0][i], sample->c[1][i], sample->c[2][i], sample->c[3][i], i);
  fclose(ccf);

  std::cerr << "sample built.\n";

  materialBase *material;
  elementBase *element;

  // UO2 TODO: Eidplacement and binding energies!
  material = new materialBase(simconf, 10.0); // rho
  element = new elementBase;
  element->z = 92; // U
  element->m = 235.0;
  element->t = 1.0;
  material->element.push_back(element);
  element = new elementBase;
  element->z = 8; // O
  element->m = 16.0;
  element->t = 2.0;
  material->element.push_back(element);
  material->prepare(); // all materials added
  sample->material.push_back(material); // add material to sample
  Real N_UO2 = material->arho;

  // xe bubble
  int gas_z1 = 54;
  material = new materialBase(simconf, 3.5); // rho
  element = new elementBase;
  element->z = gas_z1; // Xe
  element->m = 132.0;
  element->t = 1.0;
  material->element.push_back(element);
  material->prepare();
  sample->material.push_back(material); // add material to sample

  N_UO2 *= (sample->w[0]*sample->w[1]*sample->w[2] - sample->cn * 4.0/3.0 * M_PI * std::pow(r,3.0));
  std::cout << "N_UO2 = " << N_UO2 << std::endl;

  Real N_gas = sample->cn * material->arho * 4.0/3.0 * M_PI * std::pow(r,3.0);
  std::cout << "N_gas = " << N_gas << " (arho=" << material->arho << ")\n";

  // create a FIFO for recoils
  std::queue<ionBase*> recoils;

  Real norm;
  Real dif[3];

  MassInverter *m = new MassInverter;
  EnergyInverter *e = new EnergyInverter;

  Real A1, A2, Etot, E1, E2;
  int Z1, Z2;

  snprintf(fname, 199, "%s.Erec", argv[1]);
  FILE *erec = fopen(fname, "wt");

  snprintf(fname, 199, "%s.dist", argv[1]);
  FILE *rdist = fopen(fname, "wt");

  Real pos1[3];

  ionMDtag *ff1, *ff2, *pka;

  // Nev fission events
  for (int n = 0; n < Nev; n++)
  {
    if (n % 10 == 0) std::cerr << "event #" << n+1 << "\n";

    ff1 = new ionMDtag;
    ff1->gen = 0; // generation (0 = PKA)
    ff1->tag = -1;
    ff1->md = 0;

    // generate fission fragment data
    A1 = m->x(dr250());
    A2 = 235.0 - A1;
    e->setMass(A1);
    Etot = e->x(dr250());
    E1 = Etot * A2 / (A1+A2);
    E2 = Etot - E1;
    Z1 = round((A1 * 92.0) / 235.0);
    Z2 = 92 - Z1;

    ff1->_Z = Z1;
    ff1->_m = A1;
    ff1->e  = E1 * 1.0e6;

    do
    {
      for (int i = 0; i < 3; i++) ff1->dir[i] = dr250() - 0.5;
      norm = v_dot(ff1->dir, ff1->dir);
    }
    while (norm <= 0.0001);
    v_scale(ff1->dir, 1.0 / std::sqrt(norm));

    // random origin
    for (int i = 0; i < 3; i++) ff1->pos[i] = dr250() * sample->w[i];

    ff1->set_ef();
    recoils.push(ff1);

    ff2 = new ionMDtag(*ff1); // copy constructor

    // reverse direction
    for (int i = 0; i < 3; i++) ff2->dir[i] *= -1.0;

    ff2->_Z = Z2;
    ff2->_m = A2;
    ff2->e  = E2 * 1.0e6;
    ff2->md = 0;

    ff2->set_ef();
    recoils.push(ff2);

    std::cout << "A1=" << A1 << " Z1=" << Z1 << " (" << E1 << " MeV)\t"
              << "A2=" << A1 << " Z2=" << Z2 << " (" << E2 << " MeV)\n";

    // total energy of this fission event
    Real Efiss = ff1->e + ff2->e;

    while (!recoils.empty())
    {
      pka = dynamic_cast<ionMDtag*>(recoils.front());
      recoils.pop();
      sample->averages(pka);

      // do ion analysis/processing BEFORE the cascade here
      if (pka->_Z == gas_z1)
      {
	      // mark the first recoil that falls into the MD energy gap with 1
        // (child generations increase the number)
	      if (pka->e > 200 && pka->e < 12000 && pka->md == 0)
          pka->md = 1;

        if (pka->gen > 0)
        {
          // output energy and recoil generation
          fprintf(erec, "%f\t%d\t%d\n", pka->e, pka->gen, pka->md);
        }

        if (pka->tag >= 0)
        {
          for (int i = 0; i < 3; i++)
          {
            dif[i] =  sample->c[i][pka->tag] - pka->pos[i];

            if (sample->bc[i] == sampleBase::PBC)
              dif[i] -= round(dif[i] / sample->w[i]) * sample->w[i];

            pos1[i] = pka->pos[i] + dif[i];
            //printf("%f\t%f\t%f\n",   sample->c[i][pka->tag], pka->pos[i], pos1[i]);
          }
	        //printf("\n");
        }
      }

      // follow this ion's trajectory and store recoils
      // printf("%f\t%d\n", pka->e, pka->_Z);
      trim->trim(pka, recoils);

      // printf("%f %f %f %d\n", pka->pos[0], pka->pos[1], pka->pos[2], pka->tag);

      // do ion analysis/processing AFTER the cascade here

      // pka is GAS
      if (pka->_Z == gas_z1)
      {
        // output
        //printf("%f %f %f %d\n", pka->pos[0], pka->pos[1], pka->pos[2], pka->tag);

        // print out distance to cluster of origin center (and depth of recoil)
        if (pka->tag >= 0) {
          for (int i = 0; i < 3; i++)
            dif[i] = pos1[i] - pka->pos[i];

          fprintf(rdist, "%f %d %f %f %f\n", std::sqrt(v_dot(dif, dif)),
                   pka->md, pka->pos[0], pka->pos[1], pka->pos[2]);
        }
      }

      // done with this recoil
      delete pka;
    }

    // check if all energy is accounted for
    std::cout << simconf->EelTotal << std::endl;
    std::cout << simconf->EnucTotal << std::endl;
    std::cout << Efiss-(simconf->EelTotal+simconf->EnucTotal) << std::endl;
    simconf->EelTotal = 0.0;
    simconf->EnucTotal = 0.0;
  }
  fclose(rdist);
  fclose(erec);


  if (mode != PLAIN) auxout.close();

  return EXIT_SUCCESS;
}
