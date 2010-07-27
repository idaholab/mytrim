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
#include "sample_clusters.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"
#include <r250.h>

#include "functions.h"

int main(int argc, char *argv[])
{
  char fname[200];
  if( argc != 4 ) // 2
  {
    fprintf( stderr, "syntax:\nmytrim_ODS basename r Cbfactor\n\nCbfactor=1 => 1.5e-4 clusters/nm^3\n" );
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

  // initialize sample structure []
  sampleClusters *sample = new sampleClusters( 50000.0, 400.0, 400.0 );

  // initialize trim engine for the sample
  snprintf( fname, 199, "%s.phon", argv[1] );
  //FILE *phon = fopen( fname, "wt" );
  //trimPhononOut *trim = new trimPhononOut( sample, phon );
  trimBase *trim = new trimBase( sample );
  //trimBase *trim = new trimPrimaries( sample );


  //float r = 10.0;
  float r = atof( argv[2] ); //10.0;
  float Cbf = atof( argv[3] );

  sample->bc[0] = sampleBase::INF; // no PBC in x (just clusterless matrix)
  sample->initSpatialhash( int( sample->w[0] / r ) - 1,
                           int( sample->w[1] / r ) - 1,
                           int( sample->w[2] / r ) - 1 );


  // float atp = 0.1; // 10at% Mo 90at%Cu
  float v_sam = sample->w[0] * sample->w[1] * sample->w[2];
  float v_cl = 4.0/3.0 * M_PI * cub(r); 
  int n_cl; // = atp * scoef[29-1].atrho * v_sam / ( v_cl * ( ( 1.0 - atp) * scoef[42-1].atrho + atp * scoef[29-1].atrho ) );

  n_cl = v_sam * 1.5e-7 * Cbf ; // Allen08 1.5e-4/nm^3
  //fprintf( stderr, "adding %d clusters to reach %fat%% Mo\n", n_cl, atp * 100.0 );
  fprintf( stderr, "adding %d clusters...\n", n_cl );

  // cluster surfaces must be at least 25.0 Ang apart
  sample->addRandomClusters( n_cl, r, 15.0 );
  //sample->addCluster( 100.0, 100.0, 100.0, 10.0 );

  // write cluster coords with tag numbers
  snprintf( fname, 199, "%s.clcoor", argv[1] );
  FILE *ccf = fopen( fname, "wt" );
  for( int i = 0; i < sample->cn; i++)
    fprintf( ccf, "%f %f %f %f %d\n", sample->c[0][i], sample->c[1][i], sample->c[2][i], sample->c[3][i], i );
  fclose( ccf );

  fprintf( stderr, "sample built.\n" ); 
  //return 0;

  materialBase *material;
  elementBase *element;

  // Fe
  material = new materialBase( 6.98 ); // rho
  element = new elementBase;
  element->z = 26; // Fe
  element->m = 56.0;
  element->t = 1.0;
  element->Edisp = 40.0;
  material->element.push_back( element );
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
  element->z = 16; // O 
  element->m = 32.0;
  element->t = 2.0;
  material->element.push_back( element );
  material->prepare();
  sample->material.push_back( material ); // add material to sample
*/
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
  element->z = 16; // O 
  element->m = 32.0;
  element->t = 7.0;
  element->Edisp = 57.0;
  material->element.push_back( element );
  material->prepare();
  sample->material.push_back( material ); // add material to sample

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

  float A = 58.0, E = 5.0e6; int Z = 28; // 5MeV Ni

  // 100 ions
  for( int n = 0; n < 100; n++ )
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
      fprintf( erec, "%f\t%d\t%d\n", pka->e, pka->gen, pka->md );

      if( pka->z1 ==  16 || pka->z1 == 22 || pka->z1 == 39 )
      {
        if( pka->gen > 0 )
        {
          // output energy and recoil generation
          //fprintf( erec, "%f\t%d\t%d\n", pka->e, pka->gen, pka->md );
        }

        if( pka->tag >= 0 )
        {
          for( int i = 0; i < 3; i++ ) 
          {
            dif[i] =  sample->c[i][pka->tag] - pka->pos[i];
            pos2[i] = pka->pos[i];
            if( sample->bc[i] == sampleBase::PBC ) dif[i] -= round( dif[i] / sample->w[i] ) * sample->w[i];
            pos1[i] = pka->pos[i] + dif[i];
            //printf( "%f\t%f\t%f\n",   sample->c[i][pka->tag], pka->pos[i], pos1[i] );
          }
//printf( "\n" );
//if(pka->z1 == 54 && pka->gen > 0 && pka->tag >= 0 ) printf( "clust %f %f %f %d", pos1[0], pos1[1], pos1[2], pka->id );
	}
      }

      // follow this ion's trajectory and store recoils
      // printf( "%f\t%d\n", pka->e, pka->z1 );
      //pka->md = id++;

      trim->trim( pka, recoils );

      // do ion analysis/processing AFTER the cascade here

      // pka is O or Ti
      if( pka->z1 == 16 || pka->z1 == 22 || pka->z1 == 39 )
      {
        // output
        //printf( "%f %f %f %d\n", pka->pos[0], pka->pos[1], pka->pos[2], pka->tag );

        // print out distance to cluster of origin center (and depth of recoil)
        if( pka->tag >= 0 ) 
        {
          for( int i = 0; i < 3; i++ ) 
          {
            dif[i] = pos1[i] - pka->pos[i];  // distance to cluster center
            dif2[i] = pos2[i] - pka->pos[i]; // total distance it moved
          }
          fprintf( rdist, "%f %d %f %f %f %f\n", sqrt( v_dot( dif, dif ) ), pka->z1, pka->pos[0], pka->pos[1], pka->pos[2], sqrt( v_dot( dif2, dif2 ) ) );
        }


        // do a random walk
/*        jumps = 0;
        do
        {
          material = sample->lookupLayer( pka->pos );
          if( material->tag >= 0 ) break;

          do
          { 
            for( int i = 0; i < 3; i++ ) pka->dir[i] = dr250() - 0.5;
            norm = v_dot( pka->dir, pka->dir );
          }
          while( norm <= 0.0001 );
          v_scale( pka->dir, jmp / sqrtf( norm ) );

          for( int i = 0; i < 3; i++ ) pka->pos[i] += pka->dir[i];
          jumps++;
        }
        while ( pka->pos[0] > 0 && pka->pos[0] < sample->w[0] );

        if( material->tag >= 0 && jumps > 0 )
          fprintf( stderr, "walked to cluster %d (originated at %d, %d jumps)\n", material->tag, pka->tag, jumps ); */
      }

      // done with this recoil
      delete pka;

      // this should rather be done with spawnRecoil returning false
      //if( simconf->primariesOnly ) while( !recoils.empty() ) { delete recoils.front(); recoils.pop(); };
    }
  }
  fclose( rdist );
  fclose( erec );

  printf( "%d vacancies per 100 ions = %d vac/ion\n", simconf->vacancies_created, simconf->vacancies_created/100 );

  return EXIT_SUCCESS;
}