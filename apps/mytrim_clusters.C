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
#include "sample_clusters.h"
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

  // seed random number generator from system entropy pool
  FILE *urand = fopen("/dev/random", "r");
  unsigned int seed;
  if (fread(&seed, sizeof(unsigned int), 1, urand) != 1) return 1;
  fclose(urand);

  // initialize global parameter structure and read data tables from file
  SimconfType * simconf = new SimconfType(seed);
  simconf->fullTraj = false;
  simconf->tmin = 0.2;

  // initialize sample structure
  sampleClusters * sample = new sampleClusters(400.0, 400.0, 400.0);

  // initialize trim engine for the sample
  snprintf(fname, 199, "%s.phon", argv[1]);

  //FILE *phon = fopen(fname, "wt");
  //TrimPhononOut *trim = new TrimPhononOut(sample, phon);
  TrimBase *trim = new TrimBase(simconf, sample);
  //TrimBase *trim = new TrimPrimaries(sample);


  //Real r = 10.0;
  Real r = atof(argv[2]); //10.0;
  Real Cbf = atof(argv[3]);

  //sample->bc[0] = CUT; // no pbc in x dir
  sample->initSpatialhash(int(sample->w[0] / r) - 1,
                           int(sample->w[1] / r) - 1,
                           int(sample->w[2] / r) - 1);

  Real v_sam = sample->w[0] * sample->w[1] * sample->w[2];
  int n_cl;

  n_cl = v_sam * 7.0e-7 * Cbf ; // Ola06 7e-4/nm^3
  fprintf(stderr, "adding %d clusters...\n", n_cl);

  // cluster surfaces must be at least 25.0 Ang apart
  sample->addRandomClusters(n_cl, r, 15.0, simconf);

  // write cluster coords with tag numbers
  snprintf(fname, 199, "%s.clcoor", argv[1]);
  FILE *ccf = fopen(fname, "wt");
  for (int i = 0; i < sample->cn; ++i)
    fprintf(ccf, "%f %f %f %f %d\n", sample->c[0][i], sample->c[1][i], sample->c[2][i], sample->c[3][i], i);
  fclose(ccf);

  fprintf(stderr, "sample built.\n");

  MaterialBase *material;
  ElementBase *element;

  // UO2
  material = new MaterialBase(simconf, 10.0); // rho
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
  material->prepare(); // all materials added
  sample->material.push_back(material); // add material to sample

  // xe bubble
  material = new MaterialBase(simconf, 3.5); // rho
  element = new ElementBase;
  element->_Z = 54; // Xe
  element->_m = 132.0;
  element->_t = 1.0;
  // element->_t = 0.002;
  material->_element.push_back(element);
  material->prepare();
  sample->material.push_back(material); // add material to sample
  // sample->material.push_back(material); // add material to sample

  // create a FIFO for recoils
  std::queue<IonBase*> recoils;

  Real norm;
  Point dif, dif2;

  MassInverter *m = new MassInverter;
  EnergyInverter *e = new EnergyInverter;

  Real A1, A2, Etot, E1, E2;
  int Z1, Z2;

  snprintf(fname, 199, "%s.Erec", argv[1]);
  FILE *erec = fopen(fname, "wt");

  snprintf(fname, 199, "%s.dist", argv[1]);
  FILE *rdist = fopen(fname, "wt");

  Point pos1, pos2;

  IonMDTag *ff1, *ff2, *pka;

  // 1000 fission events
  for (int n = 0; n < 1000; n++) // 2000 ff
  {
    if (n % 10 == 0)
      fprintf(stderr, "pka #%d\n", n+1);

    ff1 = new IonMDTag;
    ff1->_gen = 0; // generation (0 = PKA)
    ff1->_tag = -1;
    ff1->_md = 0;
    ff1->_id = simconf->_id++;

    // generate fission fragment data
    A1 = m->x(simconf->drand());
    // A1 = 131;

    A2 = 235.0 - A1;
    Etot = e->x(simconf->drand());
    E1 = Etot * A2 / (A1 + A2);
    // E1 = 100;

    E2 = Etot - E1;
    Z1 = round((A1 * 92.0) / 235.0);
    // Z1 = 54;

    Z2 = 92 - Z1;

    ff1->_Z = Z1;
    ff1->_m = A1;
    ff1->_E  = E1 * 1.0e6;
    // ff1->_Z = 53;
    // ff1->_m = 127;
    // ff1->_E  = 70.0 * 1.0e6;

    do
    {
      for (int i = 0; i < 3; ++i) ff1->_dir(i) = simconf->drand() - 0.5;
      norm = ff1->_dir.norm_sq();
    }
    while (norm <= 0.0001 || norm > 0.25);

    /*
    norm = 1.0;
    ff1->_dir(0) = 1.0;
    ff1->_dir(1) = 0.0;
    ff1->_dir(2) = 0.0;
    */
    ff1->_dir /= std::sqrt(norm);

    // random origin (outside cluster!)
    do {
      for (int i = 0; i < 3; ++i)
        ff1->_pos(i) = simconf->drand() * sample->w[i];
    } while (sample->lookupCluster(ff1->_pos) >= 0);

    ff1->setEf();
    recoils.push(ff1);

    ff2 = new IonMDTag(*ff1); // copy constructor
    //ff1->_id = simconf->_id++;

    // reverse direction
    ff2->_dir = -ff2->_dir;

    ff2->_Z = Z2;
    ff2->_m = A2;
    ff2->_E  = E2 * 1.0e6;

    ff2->setEf();
    recoils.push(ff2);

    //fprintf(stderr, "A1=%f Z1=%d (%f MeV)\tA2=%f Z2=%d (%f MeV)\n", A1, Z1, E1, A2, Z2, E2);

    while (!recoils.empty())
    {
      pka = dynamic_cast<IonMDTag *>(recoils.front());
      recoils.pop();
      sample->averages(pka);

      // do ion analysis/processing BEFORE the cascade here

      if (pka->_Z == 54 )
      {
        // mark the first recoil that falls into the MD energy gap with 1 (child generations increase the number)
        if (pka->_E > 200 && pka->_E < 12000 && pka->_md == 0)
          pka->_md = 1;

        if (pka->_gen > 0)
        {
          // output energy and recoil generation
          fprintf(erec, "%f\t%d\t%d\n", pka->_E, pka->_gen, pka->_md);
        }

        if (pka->_tag >= 0)
        {
          for (int i = 0; i < 3; ++i)
          {
            dif(i) =  sample->c[i][pka->_tag] - pka->_pos(i);
            pos2(i) = pka->_pos(i);
            if (sample->bc[i] == SampleBase::PBC)
              dif(i) -= round(dif(i) / sample->w[i]) * sample->w[i];
            pos1(i) = pka->_pos(i) + dif(i);
            //printf("%f\t%f\t%f\n",   sample->c[i][pka->_tag], pka->_pos(i), pos1(i));
          }
          //printf("\n");
          //if (pka->_Z == 54 && pka->_gen > 0 && pka->_tag >= 0) printf("clust %f %f %f %d", pos1[0], pos1[1], pos1[2], pka->_id);
        }
      }

      // follow this ion's trajectory and store recoils
      // printf("%f\t%d\n", pka->_E, pka->_Z);
      //pka->_md = id++;

      //printf("\nstart %f %f %f %d %d %d\n", pka->_pos(0), pka->_pos(1), pka->_pos(2),  pka->_Z, pka->_md, pka->_id);
      trim->trim(pka, recoils);
      //fprintf(phon, "%f %f %f %f %d %d\n", pka->_E, pka->_pos(0), pka->_pos(1), pka->_pos(2), pka->_Z, pka->_id);

      // do ion analysis/processing AFTER the cascade here

      // pka is Xe
      if (pka->_Z == 54 )
      {
        // output
        //printf("%f %f %f %d\n", pka->_pos(0), pka->_pos(1), pka->_pos(2), pka->_tag);

        // print out distance to cluster of origin center (and depth of recoil)
        if (pka->_tag >= 0)
        {
          dif = pos1 - pka->_pos;  // distance to cluster center
          dif2 = pos2 - pka->_pos; // total distance it moved

          fprintf(rdist, "%f %d %f %f %f %f\n", dif.norm(), pka->_md, pka->_pos(0), pka->_pos(1), pka->_pos(2), dif2.norm());
        }


        // do a random walk
/*        jumps = 0;
        do
        {
          material = sample->lookupLayer(pka->_pos);
          if (material->_tag >= 0) break;

          do
          {
            for (int i = 0; i < 3; ++i) pka->_dir(i) = simconf->drand() - 0.5;
            norm = v_dot(pka->_dir, pka->_dir);
          }
          while (norm <= 0.0001);
          v_scale(pka->_dir, jmp / std::sqrt(norm));

          for (int i = 0; i < 3; ++i) pka->_pos(i) += pka->_dir(i);
          jumps++;
        }
        while (pka->_pos(0) > 0 && pka->_pos(0) < sample->w[0]);

        if (material->_tag >= 0 && jumps > 0)
          fprintf(stderr, "walked to cluster %d (originated at %d, %d jumps)\n", material->_tag, pka->_tag, jumps); */
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
