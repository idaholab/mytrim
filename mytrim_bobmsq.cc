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
#include <string.h>
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

int main(int argc, char *argv[])
{
  char fname[200];
  if( argc != 3 )
  {
    fprintf( stderr, "syntax:\n%s element energy[keV]\n",argv[0] );
    return 1;
  }

  // PKA energy
  double E = atof(argv[2])*1000.0;

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
  trimBase *trim = new trimBase( sample );
  //trimBase *trim = new trimPrimaries( sample );

  //sample->bc[0] = sampleBase::CUT; // no PBC in x (just clusterless matrix)

  // double atp = 0.1; // 10at% Mo 90at%Cu
  double v_sam = sample->w[0] * sample->w[1] * sample->w[2];
  double s_sam = sample->w[1] * sample->w[2];

  materialBase *material;
  elementBase *element;

  const char *choice[4] = {"Fe", "Si", "Cu", "Au"};
  int i;
  for(i=0; i<4 && strcmp(choice[i],argv[1])!=0; ++i);
  if(i==4) {
    fprintf( stderr, "Element choice not supported: %s\n",argv[1] );
    return 1;
  }
  double A,Z;
  switch(i) {
    case 0:
      // Fe
      material = new materialBase( 7.87 ); // rho
      element = new elementBase;
      Z = 26.0;
      A = 56.0;
      element->z = Z;
      element->m = A;
      element->t = 1.0;
      element->Edisp = 25.0;
      material->element.push_back( element );
      material->prepare(); // all materials added
      sample->material.push_back( material ); // add material to sample
      break;
    case 1:
      // Si
      material = new materialBase( 2.33 ); // rho
      element = new elementBase;
      Z = 14.0;
      A = 28.0;
      element->z = Z;
      element->m = A;
      element->t = 1.0;
      element->Edisp = 25.0;
      material->element.push_back( element );
      material->prepare(); // all materials added
      sample->material.push_back( material ); // add material to sample
      break;
    case 2:
      // Cu
      material = new materialBase( 8.89 ); // rho
      element = new elementBase;
      Z = 29.0;
      A = 63.5;
      element->z = Z;
      element->m = A;
      element->t = 1.0;
      element->Edisp = 25.0;
      material->element.push_back( element );
      material->prepare(); // all materials added
      sample->material.push_back( material ); // add material to sample
      break;
    case 3:
      // Au
      material = new materialBase( 19.32 ); // rho
      element = new elementBase;
      Z = 79.0;
      A = 197.0;
      element->z = Z;
      element->m = A;
      element->t = 1.0;
      element->Edisp = 25.0;
      material->element.push_back( element );
      material->prepare(); // all materials added
      sample->material.push_back( material ); // add material to sample
      break;
  }

  const int nstep = 1000;


  // create a FIFO for recoils
  queue<ionBase*> recoils;

  double norm;

  massInverter *m = new massInverter;
  energyInverter *e = new energyInverter;

  double A1, A2, Etot, E1, E2;
  int Z1, Z2;

  double pos1[3], pos2[3];

  ionBase *ff1, *pka;
  int id = 1;

  // squared displacement
  double sqd = 0.0, sqd2 = 0.0;

  // main loop
  for( int n = 0; n < nstep; n++ )
  {
    if( n % 10 == 0 ) fprintf( stderr, "pka #%d\n", n+1 );

    ff1 = new ionBase;
    ff1->gen = 0; // generation (0 = PKA)
    ff1->tag = -1;
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

      // store position
      if( pka->gen > 0 )
        for( int i = 0; i < 3; i++ )
          pos2[i] = pka->pos[i];

      // follow this ion's trajectory and store recoils
      trim->trim( pka, recoils );

      // do ion analysis/processing AFTER the cascade here
      if( pka->gen > 0 )
        for( int i = 0; i < 3; i++ )
          sqd += (pos2[i]-pka->pos[i])*(pos2[i]-pka->pos[i]);
      else if( pka->gen > 1 )
        for( int i = 0; i < 3; i++ )
          sqd2 += (pos2[i]-pka->pos[i])*(pos2[i]-pka->pos[i]);

      // done with this recoil
      delete pka;
    }
  }

  // output full damage data
  printf( "total sum of square displacements: %g Ang^2\n", sqd );
  printf( "%d vacancies per %d ions = %d vac/ion\n", simconf->vacancies_created, nstep, simconf->vacancies_created/nstep );
  double surf = sample->w[1] * sample->w[2];
  double natom = v_sam * sample->material[0]->arho;
  printf( "volume = %f Ang^3, surface area = %f Ang^2, containing %f atoms => %f dpa/(ion/Ang^2)\n",
          v_sam, s_sam, natom, simconf->vacancies_created / ( natom * nstep/s_sam ) );
  printf( "sqd/dpa = %g\n  sqd/vac = %g\n  sqd2/vac = %g\nnvac = %d", sqd/(simconf->vacancies_created/natom), sqd/simconf->vacancies_created, sqd2/simconf->vacancies_created,simconf->vacancies_created  );

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
