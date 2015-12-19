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
  if (argc != 4) // 2
  {
    fprintf(stderr, "syntax:\n%s basename r Cbfactor\n\nCbfactor=1 => 7e-4 bubbles/nm^3\n", argv[0]);
    return 1;
  }

  // seed randomnumber generator from system entropy pool
  FILE *urand = fopen("/dev/random", "r");
  int seed;
  fread(&seed, sizeof(int), 1, urand);
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


  //Real r = 10.0;
  Real r = atof(argv[2]); //10.0;
  Real Cbf = atof(argv[3]);


  // Real atp = 0.1; // 10at% Mo 90at%Cu
  Real v_sam = sample->w[0] * sample->w[1] * sample->w[2];
  Real v_cl = 4.0/3.0 * M_PI * cub(r);
  int n_cl; // = atp * scoef[29-1].atrho * v_sam / (v_cl * ((1.0 - atp) * scoef[42-1].atrho + atp * scoef[29-1].atrho));

  n_cl = 1;//v_sam * 7.0e-7 * Cbf ; // Ola06 7e-4/nm^3

  materialBase *material;
  elementBase *element;

  // UO2
  material = new materialBase(simconf, 9.4); // rho
  element = new elementBase;
  element->z = 92; // U
  element->m = 235.0;
  element->t = 1.0;
  material->element.push_back(element);
  element = new elementBase;
  element->z = 16; // O
  element->m = 32.0;
  element->t = 2.0;
  material->element.push_back(element);
  element = new elementBase;
  element->z = 54; // Xe
  element->m = 131.0;
  element->t = 0.0024;
  material->element.push_back(element);
  material->prepare(); // all materials added
  sample->material.push_back(material); // add material to sample

  // create a FIFO for recoils
  std::queue<ionBase*> recoils;

  Real norm;
  Real jmp = 2.7; // diffusion jump distance
  int jumps;
  Real dif[3];

  massInverter *m = new massInverter;
  energyInverter *e = new energyInverter;

  Real A1, A2, Etot, E1, E2;
  int Z1, Z2;

  snprintf(fname, 199, "%s.Erec", argv[1]);
  FILE *erec = fopen(fname, "wt");

  snprintf(fname, 199, "%s.dist", argv[1]);
  FILE *rdist = fopen(fname, "wt");

  Real pos1[3];

  ionMDtag *ff1, *ff2, *pka;
  int id = 1;

  // 5 fission events
  for (int n = 0; n < 10; n++) // 10 ff
  {
    if (n % 100 == 0) fprintf(stderr, "pka #%d\n", n+1);

    ff1 = new ionMDtag;
    ff1->gen = 0; // generation (0 = PKA)
    ff1->tag = -1;
    ff1->md = 0;
    ff1->id = simconf->id++;

    // generate fission fragment data
    A1 = m->x(dr250());
    //A1 = 131;

    A2 = 235.0 - A1;
    Etot = e->x(dr250());
    E1 = Etot * A2 / (A1 + A2);
    //E1 = 100;

    E2 = Etot - E1;
    Z1 = round((A1 * 92.0) / 235.0);
    //Z1 = 54;

    Z2 = 92 - Z1;

    /* ff1->z1 = Z1;
    ff1->m1 = A1;
    ff1->e  = E1 * 1.0e6; */

    ff1->z1 = 53;
    ff1->m1 = 127;
    ff1->e  = 70.0 * 1.0e6;

    do
    {
      for (int i = 0; i < 3; i++) ff1->dir[i] = dr250() - 0.5;
      norm = v_dot(ff1->dir, ff1->dir);
    }
    while (norm <= 0.0001);
    v_scale(ff1->dir, 1.0 / std::sqrt(norm));

    for (int i = 0; i < 3; i++) ff1->pos[i] = dr250() * sample->w[i];

    ff1->set_ef();
    recoils.push(ff1);

/*
    ff2 = new ionBase(*ff1); // copy constructor
    //ff1->id = simconf->id++;

    // reverse direction
    for (int i = 0; i < 3; i++) ff2->dir[i] *= -1.0;

    ff2->z1 = Z2;
    ff2->m1 = A2;
    ff2->e  = E2 * 1.0e6;

    ff2->set_ef();
    recoils.push(ff2);

    fprintf(stderr, "A1=%f Z1=%d (%f MeV)\tA2=%f Z2=%d (%f MeV)\n", A1, Z1, E1, A2, Z2, E2);
*/
    while (!recoils.empty())
    {
      pka = dynamic_cast<ionMDtag*>(recoils.front());
      recoils.pop();
      sample->averages(pka);

      // do ion analysis/processing BEFORE the cascade here

      if (pka->z1 == 54 )
      {
        // mark the first recoil that falls into the MD energy gap with 1 (child generations increase the number)
        if (pka->e > 200 && pka->e < 12000 && pka->md == 0) pka->md = 1;

        if (pka->gen > 0)
        {
          // output energy and recoil generation
          fprintf(erec, "%f\t%d\t%d\n", pka->e, pka->gen, pka->md);
        }

      }

      // follow this ion's trajectory and store recoils

      trim->trim(pka, recoils);

      // do ion analysis/processing AFTER the cascade here

      // pka is Xe
      if (pka->z1 == 54 )
      {
        // output
        //printf("%f %f %f %d\n", pka->pos[0], pka->pos[1], pka->pos[2], pka->tag);
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
