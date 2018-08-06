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

int
main(int argc, char * argv[])
{
  char fname[200];
  if (argc != 4) // 2
  {
    fprintf(stderr,
            "syntax:\nmytrim_ODS basename r Cbfactor\n\nCbfactor=1 => 1.5e-4 clusters/nm^3\n");
    return 1;
  }

  // seed random number generator from system entropy pool
  FILE * urand = fopen("/dev/random", "r");
  unsigned int seed;
  if (fread(&seed, sizeof(unsigned int), 1, urand) != 1)
    return 1;
  fclose(urand);

  // initialize global parameter structure and read data tables from file
  SimconfType * simconf = new SimconfType(seed);
  simconf->fullTraj = false;
  simconf->tmin = 0.2;

  // initialize sample structure []
  // sampleClusters *sample = new sampleClusters(50000.0, 400.0, 400.0);
  sampleClusters * sample = new sampleClusters(500.0, 1000.0, 1000.0);

  // initialize trim engine for the sample
  snprintf(fname, 199, "%s.phon", argv[1]);
  // FILE *phon = fopen(fname, "wt");
  // TrimPhononOut *trim = new TrimPhononOut(sample, phon);
  TrimBase * trim = new TrimBase(simconf, sample);
  // TrimBase *trim = new TrimPrimaries(sample);

  // Real r = 10.0;
  Real r = atof(argv[2]); // 10.0;
  Real Cbf = atof(argv[3]);

  sample->bc[0] = SampleBase::INF; // no PBC in x (just clusterless matrix)
  sample->initSpatialhash(
      int(sample->w[0] / r) - 1, int(sample->w[1] / r) - 1, int(sample->w[2] / r) - 1);

  // Real atp = 0.1; // 10at% Mo 90at%Cu
  Real v_sam = sample->w[0] * sample->w[1] * sample->w[2];
  // Real v_cl = 4.0/3.0 * M_PI * cub(r);
  int n_cl; // = atp * scoef[29-1].atrho * v_sam / (v_cl * ((1.0 - atp) * scoef[42-1].atrho + atp *
            // scoef[29-1].atrho));

  n_cl = v_sam * 1.5e-7 * Cbf; // Allen08 1.5e-4/nm^3
  // fprintf(stderr, "adding %d clusters to reach %fat%% Mo\n", n_cl, atp * 100.0);

  // cluster surfaces must be at least 25.0 Ang apart
  fprintf(stderr, "adding %d clusters...\n", n_cl);
  sample->addRandomClusters(n_cl, r, 15.0, simconf);

  // write cluster coords with tag numbers
  snprintf(fname, 199, "%s.clcoor", argv[1]);
  FILE * ccf = fopen(fname, "wt");
  for (int i = 0; i < sample->cn; ++i)
    fprintf(ccf,
            "%f %f %f %f %d\n",
            sample->c[0][i],
            sample->c[1][i],
            sample->c[2][i],
            sample->c[3][i],
            i);
  fclose(ccf);

  fprintf(stderr, "sample built.\n");
  // return 0;

  MaterialBase * material;
  Element element;

  /*
    // Fe
    material = new MaterialBase(7.87); // rho
    element._Z = 26; // Fe
    element._m = 56.0;
    element._t = 1.0;
    element._Edisp = 40.0;
    material->_element.push_back(element);
    material->prepare(); // all materials added
    sample->material.push_back(material); // add material to sample
  */

  // Cu
  material = new MaterialBase(simconf, 8.94); // rho
  element._Z = 29;                            // Fe
  element._m = 63.0;
  element._t = 1.0;
  element._Edisp = 40.0;
  material->_element.push_back(element);
  material->prepare();                  // all materials added
  sample->material.push_back(material); // add material to sample

  /*
    // ZrO2
    material = new MaterialBase(5.68); // rho
    element._Z = 40; // Zr
    element._m = 91.0;
    element._t = 1.0;
    material->_element.push_back(element);

    element._Z = 8; // O
    element._m = 16.0;
    element._t = 2.0;
    material->_element.push_back(element);

    material->prepare(); // all materials added
    sample->material.push_back(material); // add material to sample

    // TiO2 precipitate
    material = new MaterialBase(4.23); // rho
    element._Z = 22; // Ti
    element._m = 48.0;
    element._t = 1.0;
    material->_element.push_back(element);
    element._Z = 8; // O
    element._m = 16.0;
    element._t = 2.0;
    material->_element.push_back(element);
    material->prepare();
    sample->material.push_back(material); // add material to sample

     // Y2Ti2O7 precipitate
    material = new MaterialBase(4.6); // rho between 4.23 and 5.01
    element._Z = 39; // Y
    element._m = 89.0;
    element._t = 2.0;
    element._Edisp = 57.0;
    material->_element.push_back(element);
    element._Z = 22; // Ti
    element._m = 48.0;
    element._t = 2.0;
    element._Edisp = 57.0;
    material->_element.push_back(element);
    element._Z = 8; // O
    element._m = 16.0;
    element._t = 7.0;
    element._Edisp = 57.0;
    material->_element.push_back(element);
    material->prepare();
    sample->material.push_back(material); // add material to sample
  */
  /*
    // xe bubble
    material = new MaterialBase(3.5); // rho
    element._Z = 54; // Xe
    element._m = 132.0;
    element._t = 1.0;
    material->_element.push_back(element);
    material->prepare();
    sample->material.push_back(material); // add material to sample
  */
  // TiB2 precipitate
  material = new MaterialBase(simconf, 4.52); // rho
  element._Z = 22;                            // Ti
  element._m = 48.0;
  element._t = 1.0;
  material->_element.push_back(element);
  element._Z = 5; // B
  element._m = 11.0;
  element._t = 2.0;
  material->_element.push_back(element);
  material->prepare();
  sample->material.push_back(material); // add material to sample

  const int nstep = 1000;

  // create a FIFO for recoils
  std::queue<IonBase *> recoils;

  Real dif[3], dif2[3];

  snprintf(fname, 199, "%s.Erec", argv[1]);
  FILE * erec = fopen(fname, "wt");

  snprintf(fname, 199, "%s.dist", argv[1]);
  FILE * rdist = fopen(fname, "wt");

  Real pos1[3], pos2[3];

  IonMDTag *ff1, *pka;

  Real A = 84.0, E = 1.8e6;
  int Z = 36; // 1.8MeV Kr
  // Real A = 58.0, E = 5.0e6; int Z = 28; // 5MeV Ni
  // Real A = 56.0, E = 5.0e6; int Z = 26; // 5MeV Fe

  // 1000 ions
  for (int n = 0; n < nstep; n++)
  {
    if (n % 10 == 0)
      fprintf(stderr, "pka #%d\n", n + 1);

    ff1 = new IonMDTag;
    ff1->_gen = 0; // generation (0 = PKA)
    ff1->_tag = -1;
    ff1->_md = 0;
    ff1->_id = simconf->_id++;

    ff1->_Z = Z;
    ff1->_m = A;
    ff1->_E = E;

    ff1->_dir(0) = 1;
    ff1->_dir(1) = 0;
    ff1->_dir(2) = 0;

    ff1->_pos(0) = 0;
    ff1->_pos(1) = sample->w[1] / 2.0;
    ff1->_pos(2) = sample->w[2] / 2.0;

    ff1->setEf();
    recoils.push(ff1);

    while (!recoils.empty())
    {
      pka = dynamic_cast<IonMDTag *>(recoils.front());
      recoils.pop();
      sample->averages(pka);

      // do ion analysis/processing BEFORE the cascade here
      // fprintf(erec, "%f\t%d\t%d\n", pka->_E, pka->_gen, pka->_md);

      // pka is O or Ti
      // if (pka->_Z == 8 || pka->_Z == 22 || pka->_Z == 39)
      // pka is Xe
      // if (pka->_Z == 54)
      // pka is B or Ti
      if (pka->_Z == 5 || pka->_Z == 22)
      {
        if (pka->_gen > 0)
        {
          // output energy and recoil generation
          fprintf(erec, "%f\t%d\t%d\n", pka->_E, pka->_gen, pka->_md);
        }

        if (pka->_tag >= 0)
        {
          for (int i = 0; i < 3; ++i)
          {
            dif[i] = sample->c[i][pka->_tag] - pka->_pos(i);
            pos2[i] = pka->_pos(i);
            if (sample->bc[i] == SampleBase::PBC)
              dif[i] -= round(dif[i] / sample->w[i]) * sample->w[i];
            pos1[i] = pka->_pos(i) + dif[i];
            // printf("%f\t%f\t%f\n",   sample->c[i][pka->_tag], pka->_pos(i), pos1[i]);
          }
          // printf("\n");
          // if (pka->_Z == 54 && pka->_gen > 0 && pka->_tag >= 0) printf("clust %f %f %f %d",
          // pos1[0], pos1[1], pos1[2], pka->_id);
        }
      }

      // follow this ion's trajectory and store recoils
      // printf("%f\t%d\n", pka->_E, pka->_Z);
      // pka->_md = id++;

      trim->trim(pka, recoils);

      // do ion analysis/processing AFTER the cascade here

      // pka is O or Ti
      // if (pka->_Z == 8 || pka->_Z == 22 || pka->_Z == 39)
      // pka is Xe
      // if (pka->_Z == 54)
      // pka is B or Ti
      if (pka->_Z == 5 || pka->_Z == 22)
      {
        // output
        // printf("%f %f %f %d\n", pka->_pos(0), pka->_pos(1), pka->_pos(2), pka->_tag);

        // print out distance to cluster of origin center (and depth of recoil)
        if (pka->_tag >= 0)
        {
          for (int i = 0; i < 3; ++i)
          {
            dif[i] = pos1[i] - pka->_pos(i);  // distance to cluster center
            dif2[i] = pos2[i] - pka->_pos(i); // total distance it moved
          }
          fprintf(rdist,
                  "%f %d %f %f %f %f\n",
                  std::sqrt(v_dot(dif, dif)),
                  pka->_Z,
                  pka->_pos(0),
                  pka->_pos(1),
                  pka->_pos(2),
                  std::sqrt(v_dot(dif2, dif2)));
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
                  fprintf(stderr, "walked to cluster %d (originated at %d, %d jumps)\n",
           material->_tag, pka->_tag, jumps); */
      }

      // done with this recoil
      delete pka;

      // this should rather be done with spawnRecoil returning false
      // if (simconf->primariesOnly) while (!recoils.empty()) { delete recoils.front();
      // recoils.pop(); };
    }
  }
  fclose(rdist);
  fclose(erec);

  // output full damage data
  printf("%d vacancies per %d ions = %d vac/ion\n",
         simconf->vacancies_created,
         nstep,
         simconf->vacancies_created / nstep);
  /*
    // calculate modified kinchin pease data http://www.iue.tuwien.ac.at/phd/hoessinger/node47.html
    // just for the PKA
    Real Zatoms = 26.0, Matoms = 56.0;
    Real Epka = 5.0e6;
    Real ed = 0.0115 * std::pow(Zatoms, -7.0/3.0) * Epka;
    Real g = 3.4008 * std::pow(ed, 1.0/6.0) + 0.40244 * std::pow(ed, 3.0/4.0) + ed;
    Real kd = 0.1337 * std::pow(Zatoms, 2.0/3.0) / std::pow(Matoms, 0.5); //Z, M
    Real Ev = Epka / (1.0 + kd * g);
    Real Ed = 40.0;
    printf("%f modified PKA kinchin-pease vacancies per 100 ions = %f vac/ion\n",
            100*0.8*Ev/(2.0*Ed), 0.8*Ev/(2.0*Ed));

    // do Kinchin-Pease for all primary recoils
    printf("%f modified 1REC kinchin-pease vacancies per 100 ions = %f vac/ion\n",
            simconf->KP_vacancies, simconf->KP_vacancies / 100.0);
  */
  return EXIT_SUCCESS;
}
