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
#include <iostream>
#include <string>
#include <limits>
using namespace std;

#include "simconf.h"
#include "element.h"
#include "material.h"
#include "sample_layers.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"
#include <r250.h>

#include "functions.h"

int main(int argc, char *argv[])
{
  char fname[200];
  if( argc != 2 ) // 2
  {
    cerr << "syntax:\nmytrim_layers basename" << endl;
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
  double sx, sy, sz;
  cin >> sx >> sy >> sz;
  cin.ignore( numeric_limits<streamsize>::max(), '\n' );
  cout << "SS " << sx << ' ' << sy << ' ' << sz << endl;

  sampleLayers *sample = new sampleLayers( sx, sy, sz );
  //trimBase *trim = new trimBase( sample );
  trimBase *trim = new trimRecoils( sample );

  // Read Materials description from stdin
  int nlayer;
  double lthick, lrho, nelem;
  string lename;
  cin >> nlayer;
  cin.ignore( numeric_limits<streamsize>::max(), '\n' );
  cout << "n_layers=" << nlayer << endl;

  materialBase *material;
  elementBase *element;
  for( int i = 0; i < nlayer; i++ )
  {
    cin >> lename >> lthick >> lrho >> nelem;
    cin.ignore( numeric_limits<streamsize>::max(), '\n');
    cout << "Layer: " << lename << "  d=" << lthick << "Ang  rho=" 
         << lrho << "g/ccm  n_elements=" << nelem << endl;

    material = new materialBase( lrho ); // rho

    for( int j = 0; j < nelem; j++ )
    {
      element = new elementBase;
      cin >> lename >> element->z >> element->m >> element->t;
      cin.ignore( numeric_limits<streamsize>::max(), '\n');
      cout << "  Element: " << lename << "  Z=" << element->z 
           << "  m=" << element->m << "  fraction=" << element->t << endl;
      material->element.push_back( element );
    }

    material->prepare(); // all elements added
    sample->material.push_back( material ); // add material to sample
    sample->layerThickness.push_back( lthick );
  }

  // create a FIFO for recoils
  queue<ionBase*> recoils;

  float norm;
  float jmp = 2.7; // diffusion jump distance
  int jumps;
  float dif[3];

  massInverter *m = new massInverter;
  energyInverter *e = new energyInverter;

  float A = 74.0, E = 1.0e5;
  int Z = 36;

  snprintf( fname, 199, "%s.Erec", argv[1] );
  FILE *erec = fopen( fname, "wt" );

  snprintf( fname, 199, "%s.dist", argv[1] );
  FILE *rdist = fopen( fname, "wt" );

  float pos1[3];

  ionBase *ff1, *ff2, *pka;
  int id = 1;

  int nrec = 0;
  double sum_r2, opos[3];

  // 1000 PKA
  for( int n = 0; n < 10000000; n++ )
  {
    if( n % 1000 == 0 ) fprintf( stderr, "pka #%d\n", n+1 );

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

      //fprintf( erec, "%f\t%d\t%d\n", pka->e, pka->gen, pka->z1 );

      for( int i = 0; i < 3; i++ )
        opos[i] = pka->pos[i];

      // follow this ion's trajectory and store recoils
      trim->trim( pka, recoils );

      // do ion analysis/processing AFTER the cascade here
      if( pka->z1 != Z )
      {
        for( int i = 0; i < 3; i++ )
          sum_r2 += sqr( opos[i] - pka->pos[i] );
        nrec++;
      }

      // pka is O or Ag
      if( pka->z1 == 8 || pka->z1 == 47 ) 
      {
        // output
        //printf( "RP %f\n", pka->pos[0] );
      }

      // done with this recoil
      delete pka;

      // this should rather be done with spawnRecoil returning false
      //if( simconf->primariesOnly ) while( !recoils.empty() ) { delete recoils.front(); recoils.pop(); };
    }
  }
  fclose( rdist );
  fclose( erec );

  cout << "n=" << nrec << " sum_r2=" << sum_r2 << endl;
  return EXIT_SUCCESS;
}
