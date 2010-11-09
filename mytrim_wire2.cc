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

#include "simconf.h"
#include "element.h"
#include "material.h"
#include "sample_wire.h"
#include "sample_burried_wire.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"
#include <r250.h>

#include "functions.h"

using namespace std;

int main(int argc, char *argv[])
{
  char fname[200];
  if( argc != 8 ) 
  {
    fprintf( stderr, "syntax:\nmytrim_wire basename Eion[eV] angle[deg] numpka zpka mpka diameter(nm)\n" );
    return 1;
  }

  double epka  = atof(argv[2]);
  double theta = atof(argv[3]) * M_PI/180.0; // 0 = parallel to wire
  int numpka  = atoi(argv[4]);
  int   zpka  = atoi(argv[5]);
  double mpka  = atof(argv[6]);
  double diameter  = 10.0*atof(argv[7]);
  double length  = 10000.0;
  bool burried = true;

  // seed randomnumber generator from system entropy pool
  FILE *urand = fopen( "/dev/random", "r" );
  int seed;
  fread( &seed, sizeof(int), 1, urand );
  fclose( urand );
  r250_init( seed<0 ? -seed : seed ); // random generator goes haywire with neg. seed

  // initialize global parameter structure and read data tables from file
  simconf = new simconfType;

  // initialize sample structure
  sampleWire *sample;
  if( burried )
    sample = new sampleBurriedWire( diameter, diameter, length );
  else
  {
    sample = new sampleWire( diameter, diameter, length );
    sample->bc[2] = sampleWire::CUT;
  }
  
  // initialize trim engine for the sample
/*  const int z1 = 31;
  const int z2 = 33;
  trimVacMap *trim = new trimVacMap( sample, z1, z2 ); // GaAs*/
  trimBase *trim = new trimBase( sample );

  materialBase *material;
  elementBase *element;

  // Si
  material = new materialBase( 2.329 ); // rho
  element = new elementBase;
  element->z = 14; // Si
  element->m = 28.0;
  element->t = 1.0;
  material->element.push_back( element );
  material->prepare(); // all materials added
  sample->material.push_back( material ); // add material to sample

  // SiO2 (material[1] for the cover layer in SampleBurriedWire)
  material = new materialBase( 2.634 ); // rho
  element = new elementBase;
  element->z = 14; // Si
  element->m = 28.0;
  element->t = 1.0;
  material->element.push_back( element );
  element = new elementBase;
  element->z = 8; // O
  element->m = 16.0;
  element->t = 2.0;
  material->element.push_back( element );
  material->prepare(); // all materials added
  sample->material.push_back( material ); // add material to sample

  // create a FIFO for recoils
  queue<ionBase*> recoils;

  double norm;
  double jmp = 2.7; // diffusion jump distance
  int jumps;
  double dif[3];

  //snprintf( fname, 199, "%s.Erec", argv[1] );
  //FILE *erec = fopen( fname, "wt" );

  //snprintf( fname, 199, "%s.dist", argv[1] );
  //FILE *rdist = fopen( fname, "wt" );

  ionBase *pka;

  const int mx = 20, my = 20;
  int imap[mx][my][3];
  for( int e = 0; e < 3; e++ )
    for( int x = 0; x < mx; x++ )
      for( int y = 0; y < my; y++ )
        imap[x][y][e] = 0;

  // 10000 ions
  for( int n = 0; n < numpka; n++ )
  {
    if( n % 1000 == 0 ) fprintf( stderr, "pka #%d\n", n+1 );

    pka = new ionBase;
    pka->gen = 0; // generation (0 = PKA)
    pka->tag = -1;
    pka->md = 0;
    pka->z1 = zpka; // S
    pka->m1 = mpka;
    pka->e  = epka; 

    pka->dir[0] = 0.0;
    pka->dir[1] = sin( theta );
    pka->dir[2] = cos( theta );

    v_norm( pka->dir );

    if( burried )
    {
      // cannot anticipate the straggling in the burrial layer, thus have to shoot onto a big surface
      // TODO: take theta into account!
      pka->pos[0] = dr250() * ( 2.0*length + sample->w[0] ) - ( length + 0.5 * sample->w[0] ) ;
      pka->pos[1] = dr250() * ( 2.0*length + sample->w[1] ) - ( length + 0.5 * sample->w[1] ) ;
      pka->pos[2] = -250.0; // overcoat thickness
    }
    else
    {
      if( theta == 0.0 )
      {
        // 0 degrees => start on top of wire!
        pka->pos[2] = 0;
        do
        {
          pka->pos[0] = dr250() * sample->w[0];
          pka->pos[1] = dr250() * sample->w[1];
        } while( sample->lookupMaterial(pka->pos ) == 0 );
      }
      else
      {
        // start on side _or_ top!
        double vpos[3], t;
        do
        {
          do
          {
            vpos[0] = dr250() * sample->w[0];
            vpos[1] = 0.0;
            vpos[2] = ( dr250() * ( length + diameter/tan(theta) ) ) - diameter/tan(theta);

            t = ( 1.0 - sqrt( 1.0 - sqr( 2*vpos[0]/diameter - 1.0 ) ) ) * diameter/(2.0*pka->dir[1]);

            // if we start beyond wire length (that would be inside the substrate) then retry
          } while( t*pka->dir[2] + vpos[2] >= length );

          // if first intersection with cylinder is at z<0 then check if we hit the top face instead
          if( t*pka->dir[2] + vpos[2] < 0.0 )
            t = -vpos[2]/pka->dir[2];

          // start PKA at calculated intersection point
          for( int i = 0; i < 3; i++ )
              pka->pos[i] = t*pka->dir[i] + vpos[i];

        } while( sample->lookupMaterial(pka->pos ) == 0 );
      }
    }
    cout << "START " << pka->pos[0] << ' ' << pka->pos[1] << ' ' << pka->pos[2] << ' ' << endl;
    continue;

    pka->set_ef();
    recoils.push( pka );

    while( !recoils.empty() )
    {
      pka = recoils.front();
      recoils.pop();
      sample->averages( pka );

      // do ion analysis/processing BEFORE the cascade here

      if( pka->z1 == zpka  )
      {
        //printf(  "p1 %f\t%f\t%f\n", pka->pos[0], pka->pos[1], pka->pos[2] );
      }

      // follow this ion's trajectory and store recoils
      // printf( "%f\t%d\n", pka->e, pka->z1 );
      trim->trim( pka, recoils );

      // do ion analysis/processing AFTER the cascade here

      // ion is still in sample
      if(  sample->lookupMaterial( pka->pos ) != 0 ) 
      {
        int x, y;
        x = ( ( pka->pos[0] * mx ) / sample->w[0] );
        y = ( ( pka->pos[1] * my ) / sample->w[1] );
        x -= int(x/mx) * mx;
        y -= int(y/my) * my;

        // keep track of interstitials for the two constituents
/*        if( pka->z1 == z1 ) imap[x][y][0]++;
        else if( pka->z1 == z2 ) imap[x][y][1]++;
        else imap[x][y][2]++; // the PKAs*/
      }

      // done with this recoil
      delete pka;
    }
  }

//   char *elnam[3] = { "Ga", "As", "ion" };

//   FILE *intf, *vacf, *netf;
//   for( int e = 0; e < 3; e++ )
//   {
//     snprintf( fname, 199, "%s.%s.int", argv[1], elnam[e] );
//     intf = fopen( fname, "wt" );
//     snprintf( fname, 199, "%s.%s.vac", argv[1], elnam[e] );
//     vacf = fopen( fname, "wt" );
//     snprintf( fname, 199, "%s.%s.net", argv[1], elnam[e] );
//     netf = fopen( fname, "wt" );
// 
//     for( int y = 0; y <= my; y++ )
//     {
//       for( int x = 0; x <= mx; x++ )
//       {
//         fprintf( intf, "%d %d %d\n", x, y, (x<mx && y<my) ? imap[x][y][e] : 0 );
//         fprintf( vacf, "%d %d %d\n", x, y, (x<mx && y<my) ? trim->vmap[x][y][e] : 0 );
//         fprintf( netf, "%d %d %d\n", x, y, (x<mx && y<my) ? ( imap[x][y][e] - trim->vmap[x][y][e] ) : 0 );
//       }
//       fprintf( intf, "\n" );
//       fprintf( vacf, "\n" );
//       fprintf( netf, "\n" );
//     }
// 
//     fclose( intf );
//     fclose( vacf );
//     fclose( netf );
//   }

  return EXIT_SUCCESS;
}
