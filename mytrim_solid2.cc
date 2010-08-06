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
#include <config.h>
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
#include <r250.h>

#include "functions.h"

int main(int argc, char *argv[])
{
  char fname[200];
  if( argc != 2 )
  {
    fprintf( stderr, "syntax:\nmytrim_solid2 basename\n" );
    return 1;
  }

  // seed randomnumber generator from system entropy pool
  FILE *urand = fopen( "/dev/random", "r" );
  int seed;
  fread( &seed, sizeof(int), 1, urand );
  fclose( urand );
  r250_init( seed<0 ? -seed : seed ); // random generator goes haywire with neg. seed

  // initialize global parameter structure and read data tables from file
  simconf = new simconfType;
  simconf->fullTraj = false;
  simconf->tmin = 0.2;
  //simconf->tmin = 0.2;

  // initialize sample structure
  sampleSolid *sample = new sampleSolid( 200.0, 200.0, 200.0 );

  // initialize trim engine for the sample
  snprintf( fname, 199, "%s.phon", argv[1] );
  //FILE *phon = fopen( fname, "wt" );
  //trimPhononOut *trim = new trimPhononOut( sample, phon );
  //trimBase *trim = new trimBase( sample );
  trimBase *trim = new trimPrimaries( sample );

  sample->bc[0] = sampleBase::CUT; // no PBC in x (just clusterless matrix)

  // float atp = 0.1; // 10at% Mo 90at%Cu
  float v_sam = sample->w[0] * sample->w[1] * sample->w[2];
  float s_sam = sample->w[1] * sample->w[2];

  materialBase *material;
  elementBase *element;

/*
  // Fe
  material = new materialBase( 7.87 ); // rho
  element = new elementBase;
  element->z = 26; // Fe
  element->m = 56.0;
  element->t = 1.0;
  element->Edisp = 40.0;
  material->element.push_back( element );
  material->prepare(); // all materials added
  sample->material.push_back( material ); // add material to sample

  // ZrO2
  material = new materialBase( 5.68 ); // rho
  element = new elementBase;
  element->z = 40; // Zr
  element->m = 91.0;
  element->t = 1.0;
  material->element.push_back( element );
  element = new elementBase;
  element->z = 8; // O 
  element->m = 16.0;
  element->t = 2.0;
  material->element.push_back( element );
  material->prepare(); // all materials added
  sample->material.push_back( material ); // add material to sample
*/

  // ZrO2 Xe 0.01
  material = new materialBase( 5.68 ); // rho
  element = new elementBase;
  element->z = 40; // Zr
  element->m = 90.0;//91?
  element->t = 1.0;
  material->element.push_back( element );
  element = new elementBase;
  element->z = 8; // O 
  element->m = 16.0;
  element->t = 2.0;
  material->element.push_back( element );
/*  element = new elementBase;
  element->z = 54; // Xe 
  element->m = 132.0;
  element->t = 0.01;
  material->element.push_back( element );*/
  material->prepare(); // all materials added
  sample->material.push_back( material ); // add material to sample

/*
  // TiO2 precipitate
  material = new materialBase( 4.23 ); // rho
  element = new elementBase;
  element->z = 22; // Ti
  element->m = 48.0;
  element->t = 1.0;
  material->element.push_back( element );
  element = new elementBase;
  element->z = 8; // O 
  element->m = 16.0;
  element->t = 2.0;
  material->element.push_back( element );
  material->prepare();
  sample->material.push_back( material ); // add material to sample

   // Y2Ti2O7 precipitate
  material = new materialBase( 4.6 ); // rho between 4.23 and 5.01
  element = new elementBase;
  element->z = 39; // Y
  element->m = 89.0;
  element->t = 2.0;
  element->Edisp = 57.0;
  material->element.push_back( element );
  element = new elementBase;
  element->z = 22; // Ti
  element->m = 48.0;
  element->t = 2.0;
  element->Edisp = 57.0;
  material->element.push_back( element );
  element = new elementBase;
  element->z = 8; // O 
  element->m = 16.0;
  element->t = 7.0;
  element->Edisp = 57.0;
  material->element.push_back( element );
  material->prepare();
  sample->material.push_back( material ); // add material to sample

  // xe bubble
  material = new materialBase( 3.5 ); // rho
  element = new elementBase;
  element->z = 54; // Xe 
  element->m = 132.0;
  element->t = 1.0;
  material->element.push_back( element );
  material->prepare();
  sample->material.push_back( material ); // add material to sample
*/

  const int nstep = 10000;


  // create a FIFO for recoils
  queue<ionBase*> recoils;

  float norm;
  float jmp = 2.7; // diffusion jump distance
  int jumps;
  float dif[3], dif2[3];

  massInverter *m = new massInverter;
  energyInverter *e = new energyInverter;

  float A1, A2, Etot, E1, E2;
  int Z1, Z2;

  snprintf( fname, 199, "%s.Erec", argv[1] );
  FILE *erec = fopen( fname, "wt" );

  snprintf( fname, 199, "%s.dist", argv[1] );
  FILE *rdist = fopen( fname, "wt" );

  float pos1[3], pos2[3];

  ionBase *ff1, *pka;
  int id = 1;

  //float A = 84.0, E = 1.8e6; int Z = 36; // 1.8MeV Kr
  float A = 131.0, E = 2.0e4; int Z = 54; // 20keV Xe
  //float A = 58.0, E = 5.0e6; int Z = 28; // 5MeV Ni
  //float A = 56.0, E = 5.0e6; int Z = 26; // 5MeV Fe

  // main loop
  for( int n = 0; n < nstep; n++ )
  {
    if( n % 10 == 0 ) fprintf( stderr, "pka #%d\n", n+1 );

    ff1 = new ionBase;
    ff1->gen = 0; // generation (0 = PKA)
    ff1->tag = -1;
    ff1->md = 0;
    ff1->id = simconf->id++;

    ff1->z1 = Z;
    ff1->m1 = A;
    ff1->e  = E;

    ff1->dir[0] = 1;
    ff1->dir[1] = 0;
    ff1->dir[2] = 0;

    ff1->pos[0] = 0;
    ff1->pos[1] = sample->w[1] / 2.0;
    ff1->pos[2] = sample->w[2] / 2.0;

    ff1->set_ef();
    recoils.push( ff1 );

    while( !recoils.empty() )
    {
      pka = recoils.front();
      recoils.pop();
      sample->averages( pka );

      // do ion analysis/processing BEFORE the cascade here
      //fprintf( erec, "%f\t%d\t%d\n", pka->e, pka->gen, pka->md );

      // pka is O or Ti
      //if( pka->z1 == 8 || pka->z1 == 22 || pka->z1 == 39 )
      // pka is Xe
      float oerec = pka->e;

      if( pka->z1 == 542 )
      {
        if( pka->gen > 0 )
        {
          // output energy and recoil generation
          //fprintf( erec, "%f\t%d\t%d\n", pka->e, pka->gen, pka->md );
        }

        for( int i = 0; i < 3; i++ ) 
        {
          pos2[i] = pka->pos[i];
        }
      }

      // follow this ion's trajectory and store recoils
      // printf( "%f\t%d\n", pka->e, pka->z1 );
      //pka->md = id++;

      trim->trim( pka, recoils );
      fprintf( rdist, "%f 1\n", pka->pos[0] );

      // do ion analysis/processing AFTER the cascade here

      // pka is O or Ti
      //if( pka->z1 == 8 || pka->z1 == 22 || pka->z1 == 39 )
      // pka is Xe
      if( pka->z1 == 542 )
      {
        // output
        //printf( "%f %f %f %d\n", pka->pos[0], pka->pos[1], pka->pos[2], pka->tag );

        // print out distance to cluster of origin center (and depth of recoil)
        for( int i = 0; i < 3; i++ ) 
        {
          dif2[i] = pos2[i] - pka->pos[i]; // total distance it moved
        }
        fprintf( rdist, "%d %f %f %f %f %f\n", pka->z1, pos2[0], pos2[1], pos2[2], sqrt( v_dot( dif2, dif2 ) ), oerec );

      }

      // done with this recoil
      delete pka;

      // this should rather be done with spawnRecoil returning false
      //if( simconf->primariesOnly ) while( !recoils.empty() ) { delete recoils.front(); recoils.pop(); };
    }
  }
  fclose( rdist );
  fclose( erec );

  // output full damage data
  printf( "%d vacancies per %d ions = %d vac/ion\n", simconf->vacancies_created, nstep, simconf->vacancies_created/nstep );
  double surf = sample->w[1] * sample->w[2];
  double natom = v_sam * sample->material[0]->arho;
  printf( "volume = %f Ang^3, surface area = %f Ang^2, containing %f atoms => %f dpa/(ion/Ang^2)",
          v_sam, s_sam, natom, simconf->vacancies_created / ( natom * nstep/s_sam ) );

/*
  // calculate modified kinchin pease data http://www.iue.tuwien.ac.at/phd/hoessinger/node47.html
  // just for the PKA
  double Zatoms = 26.0, Matoms = 56.0;
  double Epka = 5.0e6;
  double ed = 0.0115 * pow( Zatoms, -7.0/3.0) * Epka;
  double g = 3.4008 * pow( ed, 1.0/6.0 ) + 0.40244 * pow( ed, 3.0/4.0 ) + ed;
  double kd = 0.1337 * pow( Zatoms, 2.0/3.0 ) / pow( Matoms, 0.5); //Z,M
  double Ev = Epka / ( 1.0 + kd * g );
  double Ed = 40.0;
  printf( "%f modified PKA kinchin-pease vacancies per 100 ions = %f vac/ion\n", 
          100*0.8*Ev/(2.0*Ed), 0.8*Ev/(2.0*Ed) );

  // do Kinchin-Pease for all primary recoils
  printf( "%f modified 1REC kinchin-pease vacancies per 100 ions = %f vac/ion\n", 
          simconf->KP_vacancies, simconf->KP_vacancies / 100.0 );
*/
  return EXIT_SUCCESS;
}
