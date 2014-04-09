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
#include <fstream>
#include <sstream>

#include "simconf.h"
#include "element.h"
#include "material.h"
#include "sample_wire.h"
#include "sample_burried_wire.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"

#include "functions.h"

using namespace std;

int main(int argc, char *argv[])
{
  char fname[200];
  if( argc != 8 )
  {
    cerr << "syntax: " << argv[0] << " basename angle[deg] diameter(nm) burried[0,1] numbermultiplier xyzout[0,1] lbinout[0,1]" << endl;
    return 1;
  }

  double theta = atof(argv[2]) * M_PI/180.0; // 0 = parallel to wire
  double diameter  = 10.0*atof(argv[3]);
  double length  = 11000.0; // 1.1 mu
  bool burried = ( atoi(argv[4]) != 0 );
  double mult = atof(argv[5]);
  bool xyz_out  = ( atoi(argv[6]) != 0 );
  bool ldat_out = ( atoi(argv[7]) != 0 );

  // ion series
  const int nstep = 5;
  double ion_dose[nstep] = { 3.0e13, 2.2e13, 1.5e13, 1.2e13, 2.5e13 }; // in ions/cm^2
  int ion_count[nstep];
  ionBase* ion_prototype[nstep];
  ion_prototype[0] = new ionBase(  5, 11.0 , 320.0e3 ); // Z,m,E
  ion_prototype[1] = new ionBase(  5, 11.0 , 220.0e3 ); // Z,m,E
  ion_prototype[2] = new ionBase(  5, 11.0 , 160.0e3 ); // Z,m,E
  ion_prototype[3] = new ionBase(  5, 11.0 , 120.0e3 ); // Z,m,E
  ion_prototype[4] = new ionBase( 15, 31.0 , 250.0e3 ); // Z,m,E

  // seed randomnumber generator from system entropy pool
  FILE *urand = fopen( "/dev/random", "r" );
  int seed;
  fread( &seed, sizeof(int), 1, urand );
  fclose( urand );
  r250_init( seed<0 ? -seed : seed ); // random generator goes haywire with neg. seed

  // initialize global parameter structure and read data tables from file
  simconf = new simconfType;
  //simconf->fullTraj = true;

  // initialize sample structure
  sampleWire *sample;
  if( burried )
    sample = new sampleBurriedWire( diameter, diameter, length );
  else
  {
    sample = new sampleWire( diameter, diameter, length );
    sample->bc[2] = sampleWire::CUT;
  }

  // calculate actual ion numbers
  for( int s = 0; s < nstep; ++s )
  {
    double A; // irradiated area in Ang^2
    if( burried )
      A =( length + sample->w[0] ) * ( length + sample->w[1] );
    else
      A = cos(theta) * M_PI * 0.25 * sample->w[0] * sample->w[1] + //   slanted top face
          sin(theta) * length * sample->w[0];                      // + projected side

    // 1cm^2 = 1e16 Ang**2, 1Ang^2 = 1e-16cm^2
    ion_count[s] = ion_dose[s] * A * 1.0e-16 * mult;
    cerr << "Ion " << s << ' ' << ion_count[s] << endl;
  }

  // initialize trim engine for the sample
/*  const int z1 = 31;
  const int z2 = 33;
  trimVacMap *trim = new trimVacMap( sample, z1, z2 ); // GaAs*/
  //trimBase *trim = new trimBase( sample );
  trimBase *trim = new trimPrimaries( sample );

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

  // map concentration along length
  int *lbins[2];
  int lx = 100; // 100 bins
  int dl = length/double(lx);
  lbins[1] = new int[lx]; // P z=15
  for( int i = 0; i < 2; ++i )
  {
    lbins[i] = new int[lx]; // 0=B (z=5), 1=P (z=15)
    for( int l = 0; l < lx; ++l )
      lbins[i][l] = 0;
  }

  // xyz data
  int xyz_lines = 0;
  stringstream xyz_data;

  for( int s = 0; s < nstep; ++s )
  {
    for( int n = 0; n < ion_count[s]; ++n )
    {
      if( n % 10000 == 0 )
        cerr << "pka #" << n+1 << endl;

      // generate new PKA from prototype ion
      pka = new ionBase( ion_prototype[s] );
      pka->gen = 0; // generation (0 = PKA)
      pka->tag = -1;

      pka->dir[0] = 0.0;
      pka->dir[1] = sin( theta );
      pka->dir[2] = cos( theta );

      v_norm( pka->dir );

      if( burried )
      {
        // cannot anticipate the straggling in the burrial layer, thus have to shoot onto a big surface
        // TODO: take theta into account!
        pka->pos[0] = ( dr250() - 0.5 ) * ( length + sample->w[0] );
        pka->pos[1] = ( dr250() - 0.5 ) * ( length + sample->w[1] );
        pka->pos[2] = -250.0; // overcoat thickness
      }
      else
      {
        if( theta == 0.0 )
        {
          // 0 degrees => start on top of wire!
          pka->pos[2] = 0.0;
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
      //cout << "START " << pka->pos[0] << ' ' << pka->pos[1] << ' ' << pka->pos[2] << ' ' << endl;
      //continue;

      pka->set_ef();
      recoils.push( pka );

      while( !recoils.empty() )
      {
        pka = recoils.front();
        recoils.pop();
        sample->averages( pka );

        // do ion analysis/processing BEFORE the cascade here

        if( pka->z1 == ion_prototype[s]->z1  )
        {
          //printf(  "p1 %f\t%f\t%f\n", pka->pos[0], pka->pos[1], pka->pos[2] );
        }

        // follow this ion's trajectory and store recoils
        trim->trim( pka, recoils );

        // do ion analysis/processing AFTER the cascade here

        // ion is in the wire
        if(  sample->lookupMaterial( pka->pos ) == sample->material[0] )
        {
          int l = pka->pos[2] / dl;
          if( l >=0 && l < lx )
          {
            if( xyz_out )
            {
              xyz_data << simconf->scoef[pka->z1-1].sym << ' '
                      << pka->pos[0]/100.0 << ' ' << pka->pos[1]/100.0 << ' ' << pka->pos[2]/100.0 << endl;
              xyz_lines++;
            }

            if( ldat_out )
              lbins[ ( pka->z1 == 5 ) ? 0 : 1 ][l]++;
          }
        }

        // done with this recoil
        delete pka;
      }
    }
  }

  // write xyz file
  if( xyz_out )
  {
    stringstream xyz_name;
    xyz_name << argv[1] << ".xyz";
    ofstream xyz( xyz_name.str().c_str() );
    xyz << xyz_lines << endl << endl << xyz_data.str();
    xyz.close();
  }

  // write lbins file (atoms per nm^3)
  if( ldat_out )
  {
    stringstream ldat_name;
    ldat_name << argv[1] << ".ldat";
    ofstream ldat( ldat_name.str().c_str() );
    double dv = 1e-3 * dl * M_PI * 0.25 *sample->w[0] * sample->w[1]; // volume per bin in nm^3
    for( int l = 0; l < lx; ++l )
      ldat << l*dl << ' ' << lbins[0][l]/(mult*dv) << ' ' << lbins[1][l]/(mult*dv) << endl;
    ldat.close();
  }
  delete[] lbins[0];
  delete[] lbins[1];

  return EXIT_SUCCESS;
}
